# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Build inputs are pinned, and the workflows run only on request with the least token access."""
from pathlib import Path
import re
import subprocess
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
import unittest
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
# GitHub runs workflows with either extension, so both are checked.
WORKFLOWS = sorted([*(ROOT / ".github/workflows").glob("*.yml"), *(ROOT / ".github/workflows").glob("*.yaml")])
# Jobs that need more than read access to the repository contents, and why.
WIDER_ACCESS = {
    ("installers.yml", "windows"): {"contents: write"},  # draft releases are visible only with write access
    ("python-publish.yml", "publish"): {"id-token: write", "contents: read"},  # trusted publishing to PyPI
    ("release.yml", "docker-cli"): {"contents: read", "actions: read"},  # checked SDK and wheel artifacts from an explicit build run
}


def jobs(text):
    """The text of each job of a workflow, by job name."""
    body = text.split("\njobs:\n", 1)[1]
    parts = re.split(r"^  ([A-Za-z0-9_-]+):\n", body, flags=re.M)
    return dict(zip(parts[1::2], parts[2::2]))


def permissions(job):
    block = re.search(r"^    permissions:\n((?:      .*\n|\s*#.*\n)+)", job, re.M)
    if block is None:
        return None
    return {line.strip() for line in block.group(1).splitlines() if line.strip() and not line.strip().startswith("#")}


class WorkflowTest(unittest.TestCase):
    def test_every_workflow_runs_only_on_request(self):
        self.assertEqual(len(WORKFLOWS), 5)
        for workflow in WORKFLOWS:
            text = workflow.read_text(encoding="utf-8")
            triggers = re.search(r"^on:\n((?:  .*\n|\n)+?)^\S", text, re.M).group(1)
            self.assertEqual(re.findall(r"^  ([a-z_]+):", triggers, re.M), ["workflow_dispatch"], workflow.name)

    def test_every_third_party_action_is_pinned_to_a_commit(self):
        for workflow in WORKFLOWS:
            for line in workflow.read_text(encoding="utf-8").splitlines():
                if "uses:" in line:
                    self.assertRegex(line, r"uses: [\w.-]+/[\w.-]+@[0-9a-f]{40} # v\d+(\.\d+)*$",
                                     workflow.name + ": " + line.strip())

    def test_each_job_states_its_token_permissions(self):
        for workflow in WORKFLOWS:
            text = workflow.read_text(encoding="utf-8")
            self.assertIn("\npermissions: {}\n", text, workflow.name)
            for name, job in jobs(text).items():
                expected = WIDER_ACCESS.get((workflow.name, name), {"contents: read"})
                self.assertEqual(permissions(job), expected, workflow.name + " " + name)

    def test_workflows_install_only_pinned_packages(self):
        # No job, and above all the one that can mint a PyPI upload token, runs a package resolved
        # at run time; the checked-out project itself bounds its extras.
        for workflow in WORKFLOWS:
            for line in workflow.read_text(encoding="utf-8").splitlines():
                if "pip install" in line:
                    for requirement in line.split("pip install", 1)[1].split():
                        if not requirement.startswith("-") and not requirement.startswith(('.', '".')):
                            self.assertRegex(requirement, r"==\d", workflow.name + ": " + line.strip())

    def test_checkouts_do_not_keep_credentials(self):
        for workflow in WORKFLOWS:
            text = workflow.read_text(encoding="utf-8")
            for checkout in re.finditer(r"uses: actions/checkout@.*\n((?:          .*\n|        with:\n)*)", text):
                self.assertIn("persist-credentials: false", checkout.group(1), workflow.name)


class PackagingTest(unittest.TestCase):
    def test_build_backend_and_optional_dependencies_have_upper_bounds(self):
        project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        requirements = project["build-system"]["requires"]
        requirements += [item for group in project["project"]["optional-dependencies"].values() for item in group]
        for requirement in requirements:
            self.assertRegex(requirement, r"<\d", requirement)

    def test_wheel_build_tools_are_pinned(self):
        config = tomllib.loads((ROOT / "scripts/cibuildwheel.toml").read_text(encoding="utf-8"))
        install = config["tool"]["cibuildwheel"]["before-build"][0]
        self.assertRegex(install, r"\bcmake==\d+\.\d+\.\d+\b")
        self.assertRegex(install, r"\bninja==\d+\.\d+\.\d+\b")

    def test_docker_base_images_are_pinned_and_the_image_is_labelled(self):
        text = (ROOT / "Dockerfile").read_text(encoding="utf-8")
        images = re.findall(r"^FROM (\S+)", text, re.M)
        self.assertEqual(len(images), 3)
        for image in images:
            self.assertRegex(image, r"^[\w.-]+:[\w.-]+@sha256:[0-9a-f]{64}$")
        version = ET.parse(ROOT / "java/pom.xml").findtext("{http://maven.apache.org/POM/4.0.0}version")
        self.assertIn('org.opencontainers.image.version="' + version + '"', text)
        self.assertIn("\nARG SMSD_REVISION\n", text)
        for label in ("title", "description", "version", "source", "licenses", "vendor"):
            self.assertIn("org.opencontainers.image." + label + "=", text)

    def test_text_files_are_checked_out_with_lf_line_endings(self):
        attributes = (ROOT / ".gitattributes").read_text(encoding="utf-8").splitlines()
        self.assertIn("* text=auto eol=lf", attributes)
        try:
            listing = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--eol"], check=True,
                                     capture_output=True, text=True).stdout
        except (OSError, subprocess.CalledProcessError):
            self.skipTest("not a git checkout")
        for line in listing.splitlines():
            index, _, attribute, path = line.split(None, 3)
            # Windows batch files are stored with CRLF and never converted.
            if index == "i/crlf":
                self.assertTrue(path.endswith((".bat", ".cmd")) and attribute == "attr/-text", line)
                continue
            self.assertIn(index, ("i/lf", "i/-text", "i/none"), line)
            if index == "i/-text":
                self.assertEqual(attribute, "attr/-text", line)

    def test_generated_wheel_and_test_caches_are_ignored(self):
        ignored = (ROOT / ".gitignore").read_text(encoding="utf-8").splitlines()
        self.assertIn("wheelhouse/", ignored)
        self.assertIn(".pytest_cache/", ignored)


if __name__ == "__main__":
    unittest.main()
