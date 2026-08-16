#!/usr/bin/env python3
"""Extract each built release bundle and prove it actually runs.

This is the gate that stands between a green build and a published release. It
deliberately tests the *archive*, not the staging tree: extraction is where the
executable bit is lost, where a CRLF shebang bites, and where a relocated
interpreter stops finding its own libraries. None of those show up in a build
log.

The bundle is extracted outside the repository and driven through its own
launcher, so nothing from the development environment can mask a missing piece.
The geometry it produces is then checked with the repository's own
tools/verify_geometry.py, run under the *development* interpreter — the bundle
is the thing under test, not the thing doing the testing.

Usage:
    python tools/smoke_bundle.py --platform win64
    python tools/smoke_bundle.py --platform linux-x86_64 --flavour portable
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INPUT_STEP = ROOT / "test" / "80mm-test-ball.step"
CC, T = 20.0, 4.0


def fail(message: str) -> int:
    print(f"  [FAIL] {message}")
    return 1


def read_version() -> str:
    """Read __version__ without importing the package."""
    text = (ROOT / "src" / "latticegen2" / "__init__.py").read_text(encoding="utf-8")
    match = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', text, re.M)
    if not match:
        raise SystemExit("FAILED: no __version__ in src/latticegen2/__init__.py")
    return match.group(1)


def find_archive(dist: Path, platform: str, flavour: str) -> Path | None:
    """Locate the archive for the version currently in the source tree.

    Matching the exact version rather than taking the newest-looking name is
    deliberate: dist/ is not cleared between builds, and sorting names as
    strings puts 2.9.0 after 2.10.0, so a "pick the last one" rule would happily
    smoke-test a stale bundle and report PASSED for something that is not the
    artifact about to be published.
    """
    version = read_version()
    for suffix in (".zip", ".tar.gz"):
        candidate = dist / f"latticegen2-{version}-{platform}-{flavour}{suffix}"
        if candidate.exists():
            return candidate
    return None


def extract(archive: Path, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    if archive.name.endswith(".zip"):
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(dest)
            # ZipFile.extractall drops permissions; restore them from the stored
            # mode so the check below tests the archive's intent, not zipfile's
            # lossiness. Real Linux users extract with tar, which keeps them.
            for info in zf.infolist():
                mode = info.external_attr >> 16
                if mode:
                    target = dest / info.filename
                    if target.exists():
                        target.chmod(mode & 0o777)
    else:
        with tarfile.open(archive, mode="r:gz") as tar:
            tar.extractall(dest)
    tops = [p for p in dest.iterdir() if p.is_dir()]
    if len(tops) != 1:
        raise SystemExit(f"FAILED: {archive.name} did not contain one top-level dir")
    return tops[0]


def check_exec_bits(bundle: Path, flavour: str) -> int:
    """On POSIX, the bundle is inert if these lost +x during archiving."""
    if os.name == "nt":
        return 0
    required = [bundle / "latticegen2.sh"]
    if flavour == "portable":
        required.append(bundle / "runtime" / "bin" / "python3")
    else:
        required.append(bundle / "install.sh")
    problems = 0
    for path in required:
        if not path.exists():
            problems += fail(f"missing {path.relative_to(bundle)}")
        elif not path.stat().st_mode & stat.S_IXUSR:
            problems += fail(f"{path.relative_to(bundle)} is not executable after extraction")
        else:
            print(f"  [ok]   {path.relative_to(bundle)} is executable")
    return problems


def check_crlf(bundle: Path) -> int:
    """Every shell script in the bundle must be LF, whatever host built it.

    Checked for all of them, not just the committed launcher: install.sh is
    generated at build time and so bypasses .gitattributes entirely.
    """
    problems = 0
    scripts = [p for p in (bundle / "latticegen2.sh", bundle / "install.sh")
               if p.exists()]
    if not any(p.name == "latticegen2.sh" for p in scripts):
        return fail("latticegen2.sh missing from bundle")
    for script in scripts:
        if b"\r\n" in script.read_bytes():
            problems += fail(
                f"{script.name} has CRLF line endings; it will not run on Linux")
        else:
            print(f"  [ok]   {script.name} has LF line endings")
    return problems


def run_install(bundle: Path) -> int:
    script = bundle / ("install.bat" if os.name == "nt" else "install.sh")
    print(f"  running {script.name}")
    result = subprocess.run(
        [str(script)] if os.name != "nt" else ["cmd", "/c", str(script)],
        cwd=bundle, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout[-3000:])
        print(result.stderr[-3000:], file=sys.stderr)
        return fail(f"{script.name} exited {result.returncode}")
    print(f"  [ok]   {script.name} installed the bundled wheels offline")
    return 0


def run_generation(bundle: Path, outdir: Path) -> tuple[int, Path | None]:
    launcher = bundle / ("latticegen2.bat" if os.name == "nt" else "latticegen2.sh")
    output = outdir / "smoke.step"
    cmd = ([str(launcher)] if os.name != "nt" else ["cmd", "/c", str(launcher)]) + [
        "-i", str(INPUT_STEP), "-cc", str(CC), "-t", str(T),
        "-o", str(output),
    ]
    print(f"  running {launcher.name} -i 80mm-test-ball.step -cc {CC} -t {T}")

    # Scrub anything that could let the bundle borrow the dev environment's
    # interpreter or packages. If the bundle only works because this machine
    # happens to have OCP installed, that is exactly what must fail here.
    env = {k: v for k, v in os.environ.items()
           if k.upper() not in {"LATTICEGEN2_PYTHON", "PYTHONPATH", "PYTHONHOME"}}

    result = subprocess.run(cmd, cwd=bundle, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        print(result.stdout[-4000:])
        print(result.stderr[-4000:], file=sys.stderr)
        return fail(f"generation exited {result.returncode}"), None

    problems = 0
    if not output.exists() or output.stat().st_size == 0:
        problems += fail("no non-empty .step was written")
    else:
        print(f"  [ok]   wrote {output.name} ({output.stat().st_size / 1e6:.1f} MB)")
    log = output.with_suffix(".log")
    if not log.exists():
        problems += fail("no .log was written next to the output")
    else:
        print(f"  [ok]   wrote {log.name}")
    return problems, (output if output.exists() else None)


def verify_geometry(step: Path) -> int:
    """Exact B-rep validity, using the repo's checks under the dev interpreter.

    A subset of the docs/specification.md §6.2 list: enough to prove the bundle
    produced real geometry rather than a plausible-looking file. The full set is
    tools/e2e.py's job and runs against the checkout, not the bundle.
    """
    sys.path.insert(0, str(ROOT / "tools"))
    import verify_geometry as vg  # noqa: E402

    problems = 0

    valid, n_solids, invalid = vg.brepcheck(str(step))
    if not valid or n_solids == 0:
        problems += fail(f"BRepCheck_Analyzer rejected {invalid} of {n_solids} solid(s)")
    else:
        print(f"  [ok]   BRepCheck_Analyzer: {n_solids} valid solid(s)")

    outside = vg.material_outside(str(step), str(INPUT_STEP))
    if abs(outside) >= (T ** 3) * 1e-6:      # same bar as tools/e2e.py
        problems += fail(f"{outside:.6g} mm^3 of material lies outside the input body")
    else:
        print(f"  [ok]   material outside input body: {outside:.6g} mm^3")

    mesh = vg.mesh_of(str(step), CC, T)
    manifold_ok, bad_edges = vg.manifold_check(mesh)
    if not manifold_ok:
        problems += fail(f"output is not a closed manifold ({bad_edges} bad edges)")
    else:
        print(f"  [ok]   closed manifold, {len(mesh.tris)} triangles")

    return problems


def smoke(platform: str, flavour: str, dist: Path, workdir: Path) -> int:
    archive = find_archive(dist, platform, flavour)
    if archive is None:
        return fail(f"no {flavour} archive for {platform} in {dist}")

    print(f"\n=== {archive.name} ({archive.stat().st_size / 1e6:.1f} MB) ===")
    bundle = extract(archive, workdir / flavour)
    print(f"  extracted to {bundle}")

    problems = check_crlf(bundle) + check_exec_bits(bundle, flavour)
    if flavour == "wheels":
        problems += run_install(bundle)

    outdir = workdir / f"out-{flavour}"
    outdir.mkdir(parents=True, exist_ok=True)
    gen_problems, step = run_generation(bundle, outdir)
    problems += gen_problems
    if step is not None:
        problems += verify_geometry(step)
    return problems


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--platform", required=True,
                        choices=["win64", "linux-x86_64"])
    parser.add_argument("--flavour", choices=["portable", "wheels", "both"],
                        default="both")
    parser.add_argument("--dist", type=Path, default=ROOT / "dist")
    parser.add_argument("--keep", action="store_true",
                        help="keep the extracted bundles for inspection")
    args = parser.parse_args(argv)

    if not INPUT_STEP.exists():
        return fail(f"smoke input {INPUT_STEP} is missing")

    flavours = ["portable", "wheels"] if args.flavour == "both" else [args.flavour]
    workdir = Path(tempfile.mkdtemp(prefix="latticegen2-smoke-"))
    print(f"smoke workspace: {workdir}")

    try:
        problems = sum(smoke(args.platform, f, args.dist, workdir) for f in flavours)
    finally:
        if args.keep:
            print(f"\nkept {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)

    print()
    if problems:
        print(f"SMOKE GATE FAILED: {problems} problem(s). Do not publish.")
        return 1
    print(f"SMOKE GATE PASSED: {', '.join(flavours)} on {args.platform}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
