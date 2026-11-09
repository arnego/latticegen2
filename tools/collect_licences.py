"""Collect the licence texts of what the pinned wheels embed, into licenses/.

Run when a pin in requirements-bundle.txt moves (docs/release.md). It needs the
network: it downloads the pinned wheels with compiled code for both release
platforms, copies the licence texts those wheels carry, and fetches from each
upstream project, at the exact tag built in, the texts a wheel leaves out.

    python tools/collect_licences.py

Two kinds of source, and the table in licenses/LICENSES.md says which each
file came from:

* **The wheel's own texts.** Pillow and NumPy append the licences of the
  libraries they embed to their wheel's licence file, per platform, so those
  files differ between Windows and Linux and both are kept.
* **Upstream, at the version compiled in.** Matplotlib's wheel carries only its
  own licence and its font licences, while it compiles in FreeType, Qhull and
  Anti-Grain Geometry and ships the fonts and assets its source tree licenses
  under ``LICENSE/``. pybind11 and cppy are header-only and compiled into
  Matplotlib, ContourPy and kiwisolver. The versions come from each project's
  build files at the pinned tag; pybind11's and cppy's exact versions are not
  recorded in a wheel, so the text is taken at the floor the build requires.

cadquery-ocp and vtk are not covered here yet: they embed many more libraries
(see licenses/LICENSES.md).
"""

from __future__ import annotations

import base64
import re
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "licenses"
GH = "https://raw.githubusercontent.com"

PLATFORMS = {
    "windows": ["--platform", "win_amd64"],
    "linux": ["--platform", "manylinux_2_28_x86_64", "--platform", "manylinux2014_x86_64",
              "--platform", "manylinux_2_17_x86_64"],
}

MATPLOTLIB_TAG = "v3.10.6"
MATPLOTLIB_ASSETS = ["AMSFONTS", "BAKOMA", "CARLOGO", "COLORBREWER", "COURIERTEN",
                     "JSXTOOLS_RESIZE_OBSERVER", "QT4_EDITOR", "SOLARIZED", "YORICK"]
UPSTREAM = {
    # file in licenses/ : URL at the tag compiled in
    "matplotlib-freetype-2.6.1-FTL.txt":
        "https://gitlab.freedesktop.org/freetype/freetype/-/raw/VER-2-6-1/docs/FTL.TXT",
    "matplotlib-qhull-8.0.2-COPYING.txt": f"{GH}/qhull/qhull/v8.0.2/COPYING.txt",
    # libavif 1.3.0 compiles these in (AVIF_CODEC_AOM / _DAV1D = LOCAL, at the
    # tags its cmake/Modules/Local*.cmake pin). Pillow's Linux wheel carries
    # their texts; its Windows wheel does not.
    "pillow-libavif-aom-3.12.1-LICENSE.txt":
        "https://aomedia.googlesource.com/aom/+/refs/tags/v3.12.1/LICENSE?format=TEXT",
    "pillow-libavif-aom-3.12.1-PATENTS.txt":
        "https://aomedia.googlesource.com/aom/+/refs/tags/v3.12.1/PATENTS?format=TEXT",
    "pillow-libavif-dav1d-1.5.1-COPYING.txt":
        "https://code.videolan.org/videolan/dav1d/-/raw/1.5.1/COPYING",
    "pybind11-LICENSE.txt": f"{GH}/pybind/pybind11/v2.13.2/LICENSE",
    "cppy-LICENSE.txt": f"{GH}/nucleic/cppy/1.3.0/LICENSE",
    **{f"matplotlib-LICENSE_{name}.txt": f"{GH}/matplotlib/matplotlib/{MATPLOTLIB_TAG}/LICENSE/LICENSE_{name}"
       for name in MATPLOTLIB_ASSETS},
}
AGG_HEADER = f"{GH}/matplotlib/matplotlib/{MATPLOTLIB_TAG}/extern/agg24-svn/include/agg_basics.h"


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as r:
        data = r.read()
    # Gitiles serves raw files base64-encoded.
    return base64.b64decode(data) if url.endswith("?format=TEXT") else data


def lf(data: bytes) -> bytes:
    return data.replace(b"\r\n", b"\n")


def pins() -> dict[str, str]:
    sys.path.insert(0, str(ROOT / "tools"))
    import build_release

    return build_release.pinned_requirements()


def download(platform: str, names: list[str], dest: Path) -> dict[str, zipfile.ZipFile]:
    p = pins()
    subprocess.run([sys.executable, "-m", "pip", "download", "-q", "--no-deps",
                    "--only-binary=:all:", "--python-version", "3.11", "--implementation", "cp",
                    *PLATFORMS[platform], "-d", str(dest),
                    *[f"{n}=={p[n]}" for n in names]], check=True)
    out = {}
    for whl in dest.glob("*.whl"):
        out[whl.name.split("-")[0].lower()] = zipfile.ZipFile(whl)
    return out


def write(name: str, data: bytes) -> None:
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(lf(data))
    print(f"  {name} ({len(data):,} bytes)")


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        for platform in PLATFORMS:
            wheels = download(platform, ["pillow", "numpy"], Path(tmp) / platform)
            write(f"pillow-LICENSE-{platform}.txt",
                  wheels["pillow"].read(next(n for n in wheels["pillow"].namelist()
                                             if n.endswith("dist-info/licenses/LICENSE"))))
            z = wheels["numpy"]
            write(f"numpy-LICENSE-{platform}.txt",
                  z.read(next(n for n in z.namelist() if n.endswith("dist-info/licenses/LICENSE.txt"))))
            if platform == "linux":
                # The licences of the code NumPy compiles in: the same on both
                # platforms but for line endings.
                for n in z.namelist():
                    m = re.search(r"dist-info/licenses/(numpy/.+[^/])$", n)
                    if m:
                        write(f"numpy-components/{m.group(1)[len('numpy/'):]}", z.read(n))
            for z in wheels.values():
                z.close()  # Windows will not remove an open file
    for name, url in UPSTREAM.items():
        write(name, fetch(url))
    header = fetch(AGG_HEADER).decode("utf-8")
    notice = header[:header.index("//----", header.index("//----") + 1)]
    write("matplotlib-agg-2.4-LICENSE.txt",
          ("Anti-Grain Geometry 2.4, as vendored by Matplotlib "
           f"{MATPLOTLIB_TAG} in extern/agg24-svn. Agg carries its licence in the "
           "header of each source file; this is the notice from "
           "include/agg_basics.h.\n\n" + notice).encode("utf-8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
