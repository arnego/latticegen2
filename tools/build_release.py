#!/usr/bin/env python3
"""Build the offline release bundles described in docs/release.md.

Two flavours are produced, for the host platform:

  portable  src + docs + licenses + a relocatable CPython with the dependencies
            already installed. Extract and run; nothing is installed on the
            target and no network is touched.

  wheels    src + docs + licenses + the dependency wheels + an install script.
            Smaller, but the target must already have a matching Python.

Windows bundles are .zip, Linux bundles are .tar.gz. That is not cosmetic: ZIP
loses the executable bit in most extractors, and the portable Linux bundle is
inert without +x on runtime/bin/python3.

Typical use::

    python tools/build_release.py --flavour both

On a machine with no network, pre-carry the inputs and point at them::

    python tools/build_release.py --flavour both \\
        --wheels-dir /media/usb/wheels --runtime /media/usb/cpython-3.11.tar.gz

Note that the payload comes from ``git archive``, so it contains committed
files only, filtered by .gitattributes ``export-ignore``. Uncommitted edits in
your working tree will NOT appear in the bundle -- commit first.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REQUIREMENTS = ROOT / "requirements-bundle.txt"

# The relocatable interpreter. The release tag is pinned so a rebuild fetches the
# same build; the patch version within it is discovered, so a new patch does not
# silently break the build. Bumping this is a deliberate edit, recorded in the
# bundle's BUNDLE-INFO.txt.
PBS_REPO = "astral-sh/python-build-standalone"
PBS_RELEASE = "20260807"
DEFAULT_TARGET_PYTHON = "3.11"

PBS_TRIPLE = {
    "win64": "x86_64-pc-windows-msvc",
    # Plain x86_64 rather than the _v2/_v3 micro-architecture variants, so the
    # bundle runs on the widest range of machines.
    "linux-x86_64": "x86_64-unknown-linux-gnu",
}

# Fixed timestamp for every archive entry, so that rebuilding a commit does not
# differ merely because time passed. 1980-01-01 is the earliest a ZIP can
# represent.
#
# This makes the wheels bundle byte-identical across rebuilds. It does NOT make
# the portable bundle so: pip generates the Windows console-script wrappers in
# runtime/Scripts/ by appending a zip to a launcher stub, and that zip carries a
# build timestamp we do not control. Measured, two clean builds of one ref: 8 of
# 7366 members differ, all of them those wrappers plus the two RECORD files
# hashing them. See docs/release.md for why they are not stripped.
ZIP_EPOCH = (1980, 1, 1, 0, 0, 0)


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def run(cmd, **kwargs):
    """Run a command, echoing it, and fail loudly on a non-zero exit."""
    printable = " ".join(str(c) for c in cmd)
    print(f"    $ {printable}", flush=True)
    return subprocess.run([str(c) for c in cmd], check=True, **kwargs)


def read_version() -> str:
    """Read __version__ without importing the package (deps may be absent)."""
    text = (ROOT / "src" / "latticegen2" / "__init__.py").read_text(encoding="utf-8")
    match = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', text, re.M)
    if not match:
        raise SystemExit("FAILED: no __version__ in src/latticegen2/__init__.py")
    return match.group(1)


def host_platform() -> str:
    if sys.platform.startswith("win"):
        return "win64"
    if sys.platform.startswith("linux"):
        return "linux-x86_64"
    raise SystemExit(
        f"FAILED: unsupported build host {sys.platform!r}. "
        "Release bundles are built for Windows and Linux x86-64 only."
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def human(size: int) -> str:
    return f"{size / 1e6:.1f} MB"


def write_text(path: Path, text: str):
    """Write a generated file with LF endings on every build host.

    Path.write_text() opens in text mode, so on Windows it silently turns every
    \\n into \\r\\n -- which would ship a CRLF install.sh. .gitattributes pins the
    committed shell scripts to LF for exactly this reason; generated ones bypass
    git entirely, so they have to say it here.
    """
    path.write_text(text, encoding="utf-8", newline="\n")


def rmtree(path: Path):
    if path.exists():
        shutil.rmtree(path, onerror=_force_remove)


def _force_remove(func, path, _exc):
    """Windows leaves read-only files (e.g. occt-LICENSE) that resist unlink."""
    os.chmod(path, stat.S_IWRITE)
    func(path)


# --------------------------------------------------------------------------
# payload staging
# --------------------------------------------------------------------------

def stage_payload(ref: str, dest: Path):
    """Extract the tracked, non-export-ignored files at `ref` into `dest`."""
    dest.mkdir(parents=True, exist_ok=True)
    archive = subprocess.run(
        ["git", "archive", "--worktree-attributes", "--format=tar", ref],
        cwd=ROOT, check=True, stdout=subprocess.PIPE,
    ).stdout
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as tar:
        tar.extractall(dest)
    for required in ("src/main.py", "src/latticegen2/gui/__init__.py"):
        if not (dest / required).exists():
            raise SystemExit(
                f"FAILED: staged payload has no {required} -- check .gitattributes")


def git_metadata(ref: str) -> tuple[str, str]:
    commit = subprocess.run(
        ["git", "rev-parse", ref], cwd=ROOT, check=True,
        stdout=subprocess.PIPE, text=True).stdout.strip()
    when = subprocess.run(
        ["git", "log", "-1", "--format=%cI", commit], cwd=ROOT, check=True,
        stdout=subprocess.PIPE, text=True).stdout.strip()
    return commit, when


# --------------------------------------------------------------------------
# dependencies
# --------------------------------------------------------------------------

def fetch_wheels(dest: Path, wheels_dir: Path | None, target_python: str) -> Path:
    """Return a directory holding the pinned wheels for this host and Python.

    A supplied --wheels-dir is used where it stands rather than copied, so
    pointing it at the default download location is harmless rather than
    self-destructive.
    """
    if wheels_dir:
        wheels_dir = wheels_dir.resolve()
        found = sorted(wheels_dir.glob("*.whl"))
        if not found:
            raise SystemExit(f"FAILED: no .whl files in {wheels_dir}")
        print(f"    using {len(found)} pre-downloaded wheels from {wheels_dir}")
        return wheels_dir

    rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    running = f"{sys.version_info.major}.{sys.version_info.minor}"
    if running != target_python:
        raise SystemExit(
            f"FAILED: wheels must be downloaded by the interpreter they target.\n"
            f"        This is Python {running}, the bundle targets {target_python}.\n"
            f"        Either run this script under Python {target_python}, pass\n"
            f"        --target-python {running}, or supply --wheels-dir."
        )
    run([sys.executable, "-m", "pip", "download", "--only-binary=:all:",
         "--disable-pip-version-check", "-r", REQUIREMENTS, "-d", dest])
    return dest


def resolve_runtime_asset(platform_name: str, target_python: str) -> tuple[str, str]:
    """Find the python-build-standalone asset for this platform in PBS_RELEASE.

    The release tag is pinned; the patch version is discovered, because pinning
    a patch that later disappears would break rebuilds of old tags for no gain.
    """
    triple = PBS_TRIPLE[platform_name]
    pattern = re.compile(
        rf"^cpython-{re.escape(target_python)}\.\d+\+{PBS_RELEASE}-"
        rf"{re.escape(triple)}-install_only_stripped\.tar\.gz$"
    )
    fallback = re.compile(
        rf"^cpython-{re.escape(target_python)}\.\d+\+{PBS_RELEASE}-"
        rf"{re.escape(triple)}-install_only\.tar\.gz$"
    )

    names: dict[str, str] = {}
    url = (f"https://api.github.com/repos/{PBS_REPO}/releases/tags/{PBS_RELEASE}")
    release = _get_json(url)
    for page in range(1, 20):
        assets = _get_json(
            f"https://api.github.com/repos/{PBS_REPO}/releases/"
            f"{release['id']}/assets?per_page=100&page={page}")
        if not assets:
            break
        for asset in assets:
            names[asset["name"]] = asset["browser_download_url"]

    for regex in (pattern, fallback):
        for name in sorted(names):
            if regex.match(name):
                return name, names[name]
    raise SystemExit(
        f"FAILED: no CPython {target_python} {triple} asset in "
        f"python-build-standalone release {PBS_RELEASE}.\n"
        f"        Update PBS_RELEASE in tools/build_release.py, or pass "
        f"--runtime with a pre-downloaded archive."
    )


def _get_json(url: str):
    request = urllib.request.Request(url, headers={
        "Accept": "application/vnd.github+json",
        "User-Agent": "latticegen2-build-release",
    })
    token = os.environ.get("GITHUB_TOKEN")
    if token:
        request.add_header("Authorization", f"Bearer {token}")
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def acquire_runtime(platform_name: str, target_python: str,
                    runtime_archive: Path | None) -> tuple[Path, str]:
    """Return a local path to the interpreter archive, and its name."""
    if runtime_archive:
        return runtime_archive, runtime_archive.name

    name, url = resolve_runtime_asset(platform_name, target_python)
    cache = ROOT / ".runtime-cache"
    cache.mkdir(exist_ok=True)
    cached = cache / name
    if cached.exists():
        print(f"    using cached runtime {cached}")
        return cached, name

    print(f"    downloading {url}")
    tmp = cached.with_suffix(".part")
    with urllib.request.urlopen(url, timeout=600) as response, tmp.open("wb") as out:
        shutil.copyfileobj(response, out)
    tmp.replace(cached)
    return cached, name


def install_runtime(stage: Path, archive: Path, wheels: Path, platform_name: str):
    """Extract the interpreter into stage/runtime and populate it from wheels."""
    runtime = stage / "runtime"
    rmtree(runtime)
    with tempfile.TemporaryDirectory() as tmpdir:
        with tarfile.open(archive, mode="r:gz") as tar:
            tar.extractall(tmpdir)
        # install_only archives contain a single top-level "python" directory.
        inner = Path(tmpdir) / "python"
        if not inner.is_dir():
            candidates = [p for p in Path(tmpdir).iterdir() if p.is_dir()]
            if len(candidates) != 1:
                raise SystemExit(
                    f"FAILED: unexpected runtime archive layout in {archive.name}")
            inner = candidates[0]
        shutil.move(str(inner), str(runtime))

    python = runtime_python(stage, platform_name)
    if not python.exists():
        raise SystemExit(f"FAILED: no interpreter at {python} after extraction")

    run([python, "-m", "pip", "install", "--no-index", "--no-warn-script-location",
         "--disable-pip-version-check",
         "--find-links", wheels, "-r", REQUIREMENTS])

    # Prove the bundled interpreter can load everything the tool imports at
    # module scope, tkinter included. `run` uses check=True, so a runtime built
    # without Tcl/Tk fails the build here rather than shipping a front-end that
    # cannot open -- and the failure names the missing module, which is the one
    # thing an operator holding a broken bundle cannot work out for themselves.
    #
    # Nothing in this repository had ever verified this: `_tkinter` comes from
    # the python-build-standalone asset, not from requirements-bundle.txt, so it
    # is outside every check the dependency pins get. On Linux it additionally
    # dlopens the system X11 libraries, which the bundle does not carry -- an
    # import here is the cheapest place to find that out.
    tk_version = subprocess.run(
        [str(python), "-c",
         "import tkinter, numpy, OCP.TopoDS; print(tkinter.TkVersion)"],
        check=True, stdout=subprocess.PIPE, text=True,
    ).stdout.strip()
    print(f"    bundled runtime loads numpy, OCP and tkinter (Tk {tk_version})")

    # Strip bytecode caches. They embed the absolute build path, so leaving them
    # in would make two builds of the same commit differ, and they regenerate on
    # the target anyway.
    removed = 0
    for cache_dir in sorted(runtime.rglob("__pycache__"), reverse=True):
        rmtree(cache_dir)
        removed += 1
    print(f"    stripped {removed} __pycache__ directories")

    collect_runtime_licenses(runtime, stage, archive.name)


def runtime_python(stage: Path, platform_name: str) -> Path:
    if platform_name == "win64":
        return stage / "runtime" / "python.exe"
    return stage / "runtime" / "bin" / "python3"


RUNTIME_LICENSE_README = """# Licenses for the bundled Python runtime

The `runtime/` directory of a portable bundle contains an unmodified CPython
distribution from python-build-standalone:

    {asset}
    https://github.com/{repo}/releases/tag/{release}

`CPYTHON-LICENSE.txt` below is that distribution's own `LICENSE.txt`, copied
verbatim. It is the PSF License Agreement together with CPython's
"Licenses and Acknowledgements for Incorporated Software" appendix, which
covers the third-party components built into the interpreter.

`pip` ships inside the runtime with its vendored dependencies' license texts
alongside their code, under
`runtime/Lib/site-packages/pip/_vendor/*/LICENSE*` (Windows) or
`runtime/lib/python*/site-packages/pip/_vendor/*/LICENSE*` (Linux). Those are
left in place rather than copied here, so each text stays next to what it
covers.

`TCL-TK-LICENSE.txt` covers Tcl/Tk, which the distribution bundles for the
standard library's `tkinter` module and which latticegen2's graphical front-end
uses. Tcl/Tk is a separate project from CPython and is not covered by the PSF
agreement or by the appendix above.

The dependencies installed into the runtime — OCCT, OCP, VTK and NumPy — are
covered by the texts in the parent directory. See `../LICENSES.md`.
"""


def collect_runtime_licenses(runtime: Path, stage: Path, asset: str):
    """Copy the interpreter distribution's own license files into licenses/runtime.

    Redistributing an interpreter means redistributing its licence documents.
    See licenses/LICENSES.md, "Bundled Python runtime".
    """
    dest = stage / "licenses" / "runtime"
    dest.mkdir(parents=True, exist_ok=True)

    # python-build-standalone's install_only(_stripped) archive puts LICENSE.txt
    # at the archive root on Windows, but under lib/python3.X/ on Linux -- the
    # root is empty there. Measured on the 3.11.15+20260807 linux-x86_64 asset:
    # no top-level LICENSE.txt, but lib/python3.11/LICENSE.txt is the same PSF
    # text. Root is checked first since that is the layout the Windows build has
    # always used.
    candidates = [runtime / "LICENSE.txt",
                  *sorted(runtime.glob("lib/python3.*/LICENSE.txt"))]
    source = next((c for c in candidates if c.is_file()), None)
    if source is not None:
        shutil.copy2(source, dest / "CPYTHON-LICENSE.txt")
        print(f"    copied runtime licence ({source.relative_to(runtime).as_posix()}) "
              f"-> licenses/runtime/CPYTHON-LICENSE.txt")
    else:
        # Not fatal to the build, but it is a compliance gap someone must close,
        # so say so loudly rather than shipping quietly.
        print("    WARNING: no LICENSE.txt in the runtime distribution; "
              "licenses/runtime/ will be incomplete -- check the archive layout")

    # Tcl/Tk is redistributed too, from the moment the front-end exists. It is
    # not a CPython component and is not covered by the PSF licence above, and
    # it is not in requirements-bundle.txt either -- it arrives inside the
    # interpreter distribution -- so it is the one dependency neither of this
    # project's two existing licence mechanisms would have caught.
    tk_licences = sorted(runtime.glob("lib/tcl*/license.terms")) + \
        sorted(runtime.glob("tcl/tcl*/license.terms")) + \
        sorted(runtime.glob("lib/tk*/license.terms"))
    if tk_licences:
        shutil.copy2(tk_licences[0], dest / "TCL-TK-LICENSE.txt")
        print(f"    copied Tcl/Tk licence ({tk_licences[0].relative_to(runtime).as_posix()}) "
              f"-> licenses/runtime/TCL-TK-LICENSE.txt")
    else:
        shutil.copy2(ROOT / "licenses" / "tcl-tk-LICENSE.txt",
                     dest / "TCL-TK-LICENSE.txt")
        print("    copied Tcl/Tk licence (this repository's own copy) "
              "-> licenses/runtime/TCL-TK-LICENSE.txt")

    write_text(dest / "README.md",
               RUNTIME_LICENSE_README.format(
                   asset=asset, repo=PBS_REPO, release=PBS_RELEASE))


# --------------------------------------------------------------------------
# generated bundle files
# --------------------------------------------------------------------------

INSTALL_BAT = """@echo off
REM Install latticegen2's dependencies from the bundled wheels. No network.
REM
REM   install.bat            create .venv beside this script (recommended)
REM   install.bat --system   install into whatever "python" resolves to
REM
REM Set LATTICEGEN2_PYTHON to choose a specific interpreter.
setlocal
set "PY=python"
if not "%LATTICEGEN2_PYTHON%"=="" set "PY=%LATTICEGEN2_PYTHON%"

REM The wheels are built for one specific Python minor version. Checking up
REM front turns an unreadable pip resolver error into a one-line explanation.
REM
REM Two separate calls, each with one unambiguous meaning. Capturing the version
REM into a variable instead would need `for /f`, which nests the command in
REM single quotes and mangles the inner quoting of "%PY%" -c "..." -- measured,
REM it silently produces no output at all. Exit codes avoid the problem, but
REM only if "cannot run" and "wrong version" are asked separately: a missing
REM executable makes cmd return 9009, which would satisfy any "errorlevel N"
REM test aimed at the version check.
"%PY%" -c "pass" >nul 2>&1
if errorlevel 1 goto :nopython
"%PY%" -c "import sys;sys.exit(0 if sys.version_info[:2]==(__MAJOR__,__MINOR__) else 3)" >nul 2>&1
if errorlevel 1 goto :badversion

if /I "%~1"=="--system" (
    "%PY%" -m pip install --no-index --find-links "%~dp0wheels" -r "%~dp0requirements-bundle.txt"
    if errorlevel 1 goto :failed
    goto :done
)

"%PY%" -m venv "%~dp0.venv"
if errorlevel 1 goto :failed
"%~dp0.venv\\Scripts\\python.exe" -m pip install --no-index --find-links "%~dp0wheels" -r "%~dp0requirements-bundle.txt"
if errorlevel 1 goto :failed

:done
echo.
echo Done. Run: latticegen2.bat -i part.step -cc 10 -t 1.5
exit /b 0

:nopython
echo FAILED: could not run a Python interpreter.
echo   Tried: %PY%
echo   Install Python __TARGET__ for x86-64 and put it on PATH, or set
echo   LATTICEGEN2_PYTHON to its full path, or use the "portable" bundle,
echo   which carries its own interpreter and needs none of this.
exit /b 1

:badversion
echo FAILED: this bundle needs Python __TARGET__.
"%PY%" -c "import sys;print('  found: Python '+sys.version.split()[0])" 2>nul
echo   The bundled wheels are built for __TARGET__ only and will not install
echo   into another version. Set LATTICEGEN2_PYTHON to a Python __TARGET__
echo   interpreter, or use the "portable" bundle, which carries its own.
exit /b 1

:failed
echo.
echo FAILED: could not install the bundled wheels from "%~dp0wheels".
exit /b 1
"""

INSTALL_SH = """#!/bin/sh
# Install latticegen2's dependencies from the bundled wheels. No network.
#
#   ./install.sh            create .venv beside this script (recommended)
#   ./install.sh --system   install into whatever python3 resolves to
#
# Set LATTICEGEN2_PYTHON to choose a specific interpreter.
set -e
DIR=$(cd "$(dirname "$0")" && pwd)
PY="${LATTICEGEN2_PYTHON:-python3}"

# The wheels are built for one specific Python minor version. Checking up front
# turns an unreadable pip resolver error into a one-line explanation.
if ! FOUND=$("$PY" -c 'import sys;print(f"{sys.version_info.major}.{sys.version_info.minor}")' 2>/dev/null); then
    echo "FAILED: could not run a Python interpreter." >&2
    echo "  Tried: $PY" >&2
    echo "  Install Python __TARGET__ for x86-64, or set LATTICEGEN2_PYTHON to its path." >&2
    exit 1
fi
if [ "$FOUND" != "__TARGET__" ]; then
    echo "FAILED: this bundle needs Python __TARGET__, but '$PY' is Python $FOUND." >&2
    echo "  The bundled wheels are built for __TARGET__ only and will not install" >&2
    echo "  into $FOUND. Set LATTICEGEN2_PYTHON to a Python __TARGET__ interpreter," >&2
    echo "  or use the 'portable' bundle instead, which carries its own." >&2
    exit 1
fi

if [ "$1" = "--system" ]; then
    "$PY" -m pip install --no-index --find-links "$DIR/wheels" -r "$DIR/requirements-bundle.txt"
else
    "$PY" -m venv "$DIR/.venv"
    "$DIR/.venv/bin/python" -m pip install --no-index \\
        --find-links "$DIR/wheels" -r "$DIR/requirements-bundle.txt"
fi

echo
echo "Done. Run: ./latticegen2.sh -i part.step -cc 10 -t 1.5"
"""


def write_install_scripts(stage: Path, target_python: str):
    major, minor = target_python.split(".")[:2]
    bat = stage / "install.bat"
    # install.bat is deliberately CRLF: cmd.exe is the consumer, and some older
    # shells mis-parse a LF-only batch file. install.sh must be LF.
    bat.write_text(
        INSTALL_BAT.replace("__TARGET__", target_python)
                   .replace("__MAJOR__", major).replace("__MINOR__", minor),
        encoding="utf-8", newline="\r\n")
    sh = stage / "install.sh"
    write_text(sh, INSTALL_SH.replace("__TARGET__", target_python))
    sh.chmod(0o755)


def write_offline_readme(stage: Path, flavour: str, platform_name: str,
                         version: str, target_python: str):
    windows = platform_name == "win64"
    launcher = "latticegen2.bat" if windows else "./latticegen2.sh"
    verify = ("certutil -hashfile <archive> SHA256"
              if windows else "sha256sum -c SHA256SUMS.txt")

    if flavour == "portable":
        setup = f"""## 1. Set up

Nothing to install. This bundle carries its own Python interpreter in
`runtime/`, with every dependency already in place.

If your extractor did not preserve permissions{'' if windows else ' (this matters on Linux)'}:

```
chmod +x latticegen2.sh runtime/bin/python3
```
"""
    else:
        setup = f"""## 1. Set up

This bundle needs **Python {target_python} on x86-64** already present. The
wheels in `wheels/` are built for exactly that version; a different minor
version will refuse to install them.

```
{'install.bat' if windows else './install.sh'}
```

That creates `.venv/` beside this file and installs from `wheels/` with
`--no-index`, so it cannot reach the network. Use `--system` instead if you
would rather install into your existing Python.
"""

    text = f"""# latticegen2 {version} — offline bundle ({flavour}, {platform_name})

Generates a parameterised diamond-strut lattice that fills the solid volume of
an input STEP file. Full documentation is in `README.md` and `docs/`.

**This bundle needs no network access at any point.**

## 0. Check it arrived intact

Compare the archive against the `SHA256SUMS.txt` published with the release:

```
{verify}
```

{setup}
## 2. Run

```
{launcher} -i part.step -cc 10 -t 1.5
```

`-cc` is the XY-plane distance between the bottom nodes of adjacent cells and
`-t` the strut profile's side length, both in mm. `-v` adds console detail; a
full `.log` is written next to the output either way. See `README.md` for the
complete parameter reference.

## 3. Things worth knowing

* **The output directory must be writable.** Boundary junctions are trimmed in
  parallel worker processes that stage intermediate `.brep` files under
  `temp/<timestamp>/` next to the output file. The folder is removed on success
  and deliberately left behind after a failure, for analysis.
* **The bundle directory itself can be read-only** once set up, including on a
  network share — output goes wherever `-o` points.
* **Every run executes at below-normal process priority**, so you can keep
  working on the machine meanwhile. Use `--cores` to hold back further.
* Exit codes are documented in `docs/algorithm.md` §10. Every failure prints one
  human-readable reason line.
{'' if flavour == 'wheels' else '''
## 4. Replacing the geometry kernel

OCCT is LGPL-2.1 and ships as stock, unmodified shared libraries inside the
`cadquery-ocp` wheel, under `runtime/`. You can replace them with a different
build of that wheel using the bundled pip, which is what preserves your rights
under that licence. See `licenses/LICENSES.md`.
'''}"""
    write_text(stage / "README-OFFLINE.md", text)


def write_bundle_info(stage: Path, *, version: str, flavour: str,
                      platform_name: str, commit: str, commit_date: str,
                      target_python: str, runtime_asset: str | None):
    lines = [
        f"latticegen2 {version}",
        f"flavour          {flavour}",
        f"platform         {platform_name}",
        f"source commit    {commit}",
        f"source date      {commit_date}",
        f"target python    {target_python}",
    ]
    if runtime_asset:
        lines.append(f"bundled runtime  {runtime_asset}")
    lines.append("")
    lines.append("Pinned dependencies (requirements-bundle.txt):")
    for line in REQUIREMENTS.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            lines.append(f"  {line}")
    lines.append("")
    write_text(stage / "BUNDLE-INFO.txt", "\n".join(lines))


# --------------------------------------------------------------------------
# deterministic archives
# --------------------------------------------------------------------------

def iter_files(root: Path):
    """Yield (absolute path, archive-relative posix path), sorted, dirs first."""
    entries = []
    for path in root.rglob("*"):
        entries.append((path, path.relative_to(root).as_posix()))
    entries.sort(key=lambda item: item[1])
    return entries


def is_executable(path: Path, rel: str) -> bool:
    if rel.endswith(".sh"):
        return True
    if rel.startswith("runtime/bin/"):
        return True
    try:
        return bool(path.stat().st_mode & stat.S_IXUSR)
    except OSError:
        return False


def make_zip(stage: Path, top: str, out: Path):
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED,
                         compresslevel=6) as zf:
        for path, rel in iter_files(stage):
            arcname = f"{top}/{rel}"
            if path.is_dir():
                info = zipfile.ZipInfo(arcname + "/", date_time=ZIP_EPOCH)
                info.external_attr = (0o40755 << 16) | 0x10
                zf.writestr(info, b"")
                continue
            info = zipfile.ZipInfo(arcname, date_time=ZIP_EPOCH)
            info.compress_type = zipfile.ZIP_DEFLATED
            mode = 0o755 if is_executable(path, rel) else 0o644
            info.external_attr = mode << 16
            # Streamed rather than read whole: the bundle contains an 87.6 MB
            # OCP extension module and comparable VTK libraries, and writestr
            # would hold each one entirely in memory.
            with path.open("rb") as src, zf.open(info, "w") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)


def make_targz(stage: Path, top: str, out: Path):
    with out.open("wb") as raw:
        # mtime=0 keeps the gzip header itself free of a build timestamp.
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            with tarfile.open(fileobj=gz, mode="w") as tar:
                for path, rel in iter_files(stage):
                    info = tar.gettarinfo(str(path), arcname=f"{top}/{rel}")
                    info.mtime = 0
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    if info.isdir():
                        info.mode = 0o755
                        tar.addfile(info)
                    elif info.issym():
                        tar.addfile(info)
                    else:
                        info.mode = 0o755 if is_executable(path, rel) else 0o644
                        with path.open("rb") as handle:
                            tar.addfile(info, handle)


# --------------------------------------------------------------------------
# orchestration
# --------------------------------------------------------------------------

def build_one(flavour: str, *, version: str, platform_name: str, ref: str,
              out_dir: Path, runtime_archive: Path | None,
              target_python: str, shared_wheels: Path) -> Path:
    top = f"latticegen2-{version}-{platform_name}-{flavour}"
    stage = out_dir / "stage" / top
    rmtree(stage)

    print(f"[{flavour}] staging payload from {ref}")
    stage_payload(ref, stage)
    commit, commit_date = git_metadata(ref)

    runtime_asset = None
    if flavour == "portable":
        print(f"[{flavour}] acquiring relocatable CPython {target_python}")
        archive, runtime_asset = acquire_runtime(
            platform_name, target_python, runtime_archive)
        print(f"[{flavour}] installing dependencies into runtime/")
        install_runtime(stage, archive, shared_wheels, platform_name)
        # The wheels are deliberately NOT duplicated here: the same content is
        # already installed under runtime/, and carrying both would add
        # 140-230 MB to the asset for no benefit. The wheels flavour, published
        # alongside, is where they live.
    else:
        print(f"[{flavour}] copying wheels")
        dest = stage / "wheels"
        dest.mkdir(parents=True, exist_ok=True)
        for wheel in sorted(shared_wheels.glob("*.whl")):
            shutil.copy2(wheel, dest / wheel.name)
        write_install_scripts(stage, target_python)

    write_offline_readme(stage, flavour, platform_name, version, target_python)
    write_bundle_info(stage, version=version, flavour=flavour,
                      platform_name=platform_name, commit=commit,
                      commit_date=commit_date, target_python=target_python,
                      runtime_asset=runtime_asset)

    out_dir.mkdir(parents=True, exist_ok=True)
    if platform_name == "win64":
        out = out_dir / f"{top}.zip"
        out.unlink(missing_ok=True)
        make_zip(stage, top, out)
    else:
        out = out_dir / f"{top}.tar.gz"
        out.unlink(missing_ok=True)
        make_targz(stage, top, out)

    print(f"[{flavour}] {out.name}  {human(out.stat().st_size)}")
    print(f"[{flavour}] sha256 {sha256_file(out)}")
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Build latticegen2 offline release bundles.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Typical use::", 1)[1] if "Typical use::" in __doc__ else None,
    )
    parser.add_argument("--flavour", choices=["portable", "wheels", "both"],
                        default="both", help="which bundle(s) to build (default: both)")
    parser.add_argument("--ref", default="HEAD",
                        help="git ref to package (default: HEAD). Uncommitted "
                             "working-tree changes are never included.")
    parser.add_argument("--out", type=Path, default=ROOT / "dist",
                        help="output directory (default: dist/)")
    parser.add_argument("--wheels-dir", type=Path, default=None,
                        help="use pre-downloaded wheels instead of fetching them")
    parser.add_argument("--runtime", type=Path, default=None,
                        help="use a pre-downloaded python-build-standalone archive")
    parser.add_argument("--target-python", default=DEFAULT_TARGET_PYTHON,
                        help=f"Python minor version the bundle targets "
                             f"(default: {DEFAULT_TARGET_PYTHON})")
    parser.add_argument("--keep-stage", action="store_true",
                        help="leave the staging tree in place for inspection")
    args = parser.parse_args(argv)

    version = read_version()
    platform_name = host_platform()
    flavours = ["portable", "wheels"] if args.flavour == "both" else [args.flavour]

    print(f"latticegen2 {version} -> {platform_name}, flavours: {', '.join(flavours)}")

    # One wheel set shared by both flavours.
    print("[deps] resolving wheels")
    shared_wheels = fetch_wheels(args.out / "wheels", args.wheels_dir,
                                 args.target_python)
    wheels = sorted(shared_wheels.glob("*.whl"))
    total = sum(w.stat().st_size for w in wheels)
    print(f"[deps] {len(wheels)} wheels, {human(total)}")

    built = [
        build_one(flavour, version=version, platform_name=platform_name,
                  ref=args.ref, out_dir=args.out,
                  runtime_archive=args.runtime, target_python=args.target_python,
                  shared_wheels=shared_wheels)
        for flavour in flavours
    ]

    if not args.keep_stage:
        rmtree(args.out / "stage")

    print("\nBuilt:")
    for path in built:
        print(f"  {path.name:60s} {human(path.stat().st_size):>10s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
