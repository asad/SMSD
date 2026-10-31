# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Copy the licence files of the installed smsd wheel to a documentation folder.

The wheel carries LICENSE, NOTICE and the third-party licences of the code it
contains; the image keeps copies under /usr/share/doc/smsd, as the CPU image does.
"""
from importlib.metadata import distribution
from pathlib import Path, PurePosixPath
import sys

REQUIRED = ("LICENSE", "NOTICE")


def licence_files(dist):
    """Pairs of (installed file, path below the documentation folder)."""
    files = [PurePosixPath(str(file)) for file in dist.files or ()]
    info = [file for file in files if len(file.parts) > 1 and file.parts[0].endswith(".dist-info")]
    # Wheels following PEP 639 keep licence files below <name>.dist-info/licenses/.
    nested = any(file.parts[1:] == ("licenses", "LICENSE") for file in info)
    for file in info:
        if nested:
            if file.parts[1] != "licenses":
                continue
            inner = file.parts[2:]
        else:
            inner = file.parts[1:]
        if inner and (inner[0] in REQUIRED and len(inner) == 1
                      or inner[0] == "licenses" and len(inner) > 1):
            yield file, PurePosixPath(*inner)


def copy_licences(destination, dist=None):
    """Copy the files and return their paths relative to destination."""
    dist = dist or distribution("smsd")
    copied = []
    for file, relative in licence_files(dist):
        target = Path(destination, relative)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(Path(dist.locate_file(file)).read_bytes())
        copied.append(relative.as_posix())
    missing = [name for name in REQUIRED if name not in copied]
    if missing:
        raise ValueError("The installed smsd wheel has no " + " or ".join(missing) + " file")
    return sorted(copied)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        raise SystemExit("usage: copy_licences.py DESTINATION")
    for name in copy_licences(argv[0]):
        print("Copied", name)


if __name__ == "__main__":
    main()
