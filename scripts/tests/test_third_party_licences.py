# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Third-party licence files and SPDX metadata in the CUDA, conda and container recipes."""
import hashlib
import importlib.util
from importlib.metadata import PathDistribution
from pathlib import Path
import re
import tempfile
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[2]
CUDA = ROOT / "packaging/smsd-cuda"
# The CCCL top-level LICENSE of NVIDIA/cccl release 2.7.0 (commit b5fe509f), byte for byte.
CCCL_SIZE = 33838
CCCL_SHA256 = "01b767dcd7d36f42efb608076741cf83f154a995e198028cb698aadc3a43b63b"
CUDA_LICENCE = "Apache-2.0 AND BSD-3-Clause AND Apache-2.0 WITH LLVM-exception"
LICENCES = {"Apache-2.0", "BSD-3-Clause", "PSF-2.0", "LicenseRef-NVIDIA-CUDA"}
EXCEPTIONS = {"LLVM-exception"}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def valid_expression(expression):
    """A flat SPDX AND expression of known identifiers, each optionally WITH an exception."""
    terms = expression.split(" AND ")
    for term in terms:
        match = re.fullmatch(r"([A-Za-z0-9.+-]+)(?: WITH ([A-Za-z0-9.+-]+))?", term)
        if not match or match.group(1) not in LICENCES or (
                match.group(2) is not None and match.group(2) not in EXCEPTIONS):
            return False
    return len(set(terms)) == len(terms)


class CcclLicence(unittest.TestCase):
    def test_cccl_licence_is_the_upstream_file(self):
        data = (CUDA / "licenses/cccl/LICENSE").read_bytes()
        self.assertEqual(len(data), CCCL_SIZE)
        self.assertEqual(hashlib.sha256(data).hexdigest(), CCCL_SHA256)
        text = data.decode("utf-8")
        self.assertIn("libcu++ is under the Apache License v2.0 with LLVM Exceptions", text)
        self.assertIn("CUB's source code is released under the BSD 3-Clause license", text)

    def test_cccl_licence_is_installed_and_listed_beside_cub(self):
        cmake = (CUDA / "CMakeLists.txt").read_text()
        for name in ("cub", "cccl"):
            self.assertIn(f'install(FILES "${{CMAKE_CURRENT_SOURCE_DIR}}/licenses/{name}/LICENSE"\n'
                          f'    DESTINATION "${{CMAKE_INSTALL_DATADIR}}/licenses/smsd-cuda/licenses/{name}")',
                          cmake)
        recipe = (ROOT / "packaging/conda/cuda/meta.yaml").read_text()
        for name in ("cub", "cccl"):
            self.assertIn(f"    - packaging/smsd-cuda/licenses/{name}/LICENSE\n", recipe)
        dockerfile = (ROOT / "Dockerfile.cuda").read_text()
        ignored = (ROOT / ".dockerignore").read_text().splitlines()
        for name in ("cub", "cccl"):
            self.assertIn(f"COPY licenses/{name}/LICENSE ./licenses/{name}/LICENSE", dockerfile)
            self.assertIn(f"/usr/share/doc/smsd/licenses/{name}/LICENSE", dockerfile)
            # The Docker context excludes licenses/* apart from the files named after it.
            self.assertGreater(ignored.index(f"!licenses/{name}/LICENSE"), ignored.index("licenses/*"))
        self.assertIn("licenses/cccl/LICENSE", (ROOT / "scripts/docker_install.py").read_text())
        self.assertIn("licenses/cccl/LICENSE", (CUDA / "NOTICE").read_text())


class LicenceMetadata(unittest.TestCase):
    def test_cuda_metadata_names_the_cub_and_cccl_terms(self):
        project = tomllib.loads((CUDA / "pyproject.toml").read_text())["project"]
        self.assertEqual(project["license"], {"text": CUDA_LICENCE})
        recipe = (ROOT / "packaging/conda/cuda/meta.yaml").read_text()
        self.assertEqual(re.search(r"^  license: (.+)$", recipe, re.M).group(1), CUDA_LICENCE)
        label = re.search(r'org\.opencontainers\.image\.licenses="([^"]+)"',
                          (ROOT / "Dockerfile.cuda").read_text()).group(1)
        self.assertTrue(label.startswith(CUDA_LICENCE + " AND "), label)
        for expression in (CUDA_LICENCE, label):
            self.assertTrue(valid_expression(expression), expression)

    def test_main_conda_recipe_lists_the_licences_of_bundled_code(self):
        recipe = (ROOT / "packaging/conda/meta.yaml").read_text()
        section = recipe.split("  license_file:\n", 1)[1].split("\n  summary:", 1)[0]
        listed = re.findall(r"^    - (\S+)$", section, re.M)
        wheel = tomllib.loads((ROOT / "pyproject.toml").read_text())["tool"]["scikit-build"]["wheel"]["license-files"]
        # The conda build links no OpenMP runtime and bundles no MSVC runtime.
        expected = [name for name in wheel if not name.startswith(("licenses/libomp/", "licenses/libgomp/",
                                                                    "licenses/msvc/"))]
        self.assertEqual(listed, expected)
        for name in listed:
            self.assertTrue((ROOT / name).is_file(), name)


class McpImageLicences(unittest.TestCase):
    COPY = load("mcp_copy_licences", ROOT / "packaging/mcp/copy_licences.py")

    def distribution(self, root, nested):
        info = Path(root) / "site/smsd-8.1.0.dist-info"
        names = ["LICENSE", "NOTICE", "licenses/rdkit/LICENSE", "licenses/pybind11/LICENSE"]
        names = ["licenses/" + name for name in names] if nested else names
        for name in names + ["METADATA"]:
            (info / name).parent.mkdir(parents=True, exist_ok=True)
            (info / name).write_text(name + "\n")
        (info / "RECORD").write_text("".join(f"smsd-8.1.0.dist-info/{name},,\n" for name in names + ["METADATA"])
                                     + "smsd/__init__.py,,\n")
        return PathDistribution(info)

    def test_wheel_licences_are_copied_for_both_dist_info_layouts(self):
        for nested in (True, False):
            with self.subTest(nested=nested), tempfile.TemporaryDirectory() as root:
                destination = Path(root) / "doc"
                copied = self.COPY.copy_licences(destination, self.distribution(root, nested))
                self.assertEqual(copied, ["LICENSE", "NOTICE", "licenses/pybind11/LICENSE",
                                          "licenses/rdkit/LICENSE"])
                prefix = "licenses/" if nested else ""
                self.assertEqual((destination / "licenses/rdkit/LICENSE").read_text(),
                                 prefix + "licenses/rdkit/LICENSE\n")

    def test_missing_notice_is_an_error(self):
        with tempfile.TemporaryDirectory() as root:
            dist = self.distribution(root, True)
            record = Path(root) / "site/smsd-8.1.0.dist-info/RECORD"
            record.write_text(record.read_text().replace("smsd-8.1.0.dist-info/licenses/NOTICE,,\n", ""))
            with self.assertRaisesRegex(ValueError, "NOTICE"):
                self.COPY.copy_licences(Path(root) / "doc", dist)

    def test_image_runs_the_copy_after_installing_the_wheel(self):
        dockerfile = (ROOT / "packaging/mcp/Dockerfile").read_text()
        self.assertIn("COPY copy_licences.py /opt/copy_licences.py", dockerfile)
        install = dockerfile.index("python -m pip install")
        self.assertGreater(dockerfile.index("python /opt/copy_licences.py /usr/share/doc/smsd"), install)
        # The MCP build context excludes everything apart from the files named after "*".
        ignored = (ROOT / "packaging/mcp/.dockerignore").read_text().splitlines()
        self.assertGreater(ignored.index("!copy_licences.py"), ignored.index("*"))


if __name__ == "__main__":
    unittest.main()
