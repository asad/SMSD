# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Release preparation must preserve an existing collection of assets."""
import os
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import venv


@unittest.skipUnless(os.name == "posix" and shutil.which("bash"), "Unix release preparation script")
class PrepareReleaseTest(unittest.TestCase):
    def test_existing_output_is_rejected_before_building(self):
        script = Path(__file__).resolve().parents[1] / "prepare-release.sh"
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            runtime = directory / "runtime"
            venv.EnvBuilder(with_pip=False, symlinks=True).create(runtime)
            commands = directory / "commands"
            commands.mkdir()
            maven = commands / "mvn"
            maven.write_text("#!/bin/sh\nprintf '%s\\n' 'Java version: 25.0.3'\n")
            maven.chmod(0o755)
            output = directory / "release"
            output.mkdir()
            sentinel = output / "collected-installer.msi"
            sentinel.write_bytes(b"previously validated installer")
            environment = dict(os.environ, SMSD_RELEASE_PYTHON=str(runtime / "bin/python"),
                               SMSD_RELEASE_OUTPUT_DIR=str(output), PATH=str(commands) + os.pathsep + os.environ["PATH"])
            result = subprocess.run(["bash", str(script)], env=environment,
                                    text=True, capture_output=True, timeout=20)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Release output already exists", result.stderr)
            self.assertEqual(sentinel.read_bytes(), b"previously validated installer")
            self.assertEqual(list(output.iterdir()), [sentinel])

    def test_failed_preparation_retains_the_canonical_source_and_staging(self):
        script = Path(__file__).resolve().parents[1] / "prepare-release.sh"
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            checkout = directory / "checkout"
            (checkout / "scripts").mkdir(parents=True)
            shutil.copy2(script, checkout / "scripts/prepare-release.sh")
            commands = directory / "commands"
            commands.mkdir()
            calls = directory / "calls.jsonl"
            python = commands / "python"
            python.write_text("#!" + sys.executable + "\n" + '''
import json, os, sys
from pathlib import Path
with Path(os.environ["CALLS"]).open("a") as stream:
    stream.write(json.dumps(sys.argv[1:]) + "\\n")
if sys.argv[1] == "-c" and "ElementTree" in sys.argv[2]:
    print("8.0.3")
if sys.argv[1] == "-c" and "sys.stdin" in sys.argv[2]:
    sys.stdin.read()
if sys.argv[1:3] == ["-m", "pytest"]:
    sys.exit(1)
''')
            python.chmod(0o755)
            maven = commands / "mvn"
            maven.write_text("#!/bin/sh\nprintf '%s\\n' 'Java version: 25.0.3'\n")
            maven.chmod(0o755)
            java = commands / "java"
            java.write_text('#!/bin/sh\nprintf \'%s\\n\' \'java version "1.8.0_504"\' >&2\n')
            java.chmod(0o755)
            (commands / "bin").mkdir()
            (commands / "bin/java").symlink_to(java)
            canonical = directory / "canonical"
            canonical.mkdir()
            source = canonical / "smsd-8.0.3.tar.gz"
            source.write_bytes(b"unchanged canonical source")
            manifest = canonical / "SOURCE_MANIFEST.json"
            manifest.write_bytes(b"unchanged canonical manifest")
            environment = dict(os.environ, SMSD_RELEASE_PYTHON=str(python), CALLS=str(calls),
                               SMSD_RELEASE_SDIST=str(source), SMSD_RELEASE_COMMIT="a" * 40,
                               SMSD_RELEASE_JAVA8_HOME=str(commands),
                               PATH=str(commands) + os.pathsep + os.environ["PATH"])
            result = subprocess.run(["bash", str(checkout / "scripts/prepare-release.sh")],
                                    env=environment, text=True, capture_output=True, timeout=20)
            self.assertEqual(result.returncode, 1, result.stderr)
            assets, = (checkout / "build").glob("release-assets.*")
            source_directory, = (checkout / "build").glob("release-source.*")
            self.assertTrue(source_directory.is_dir())
            self.assertEqual((assets / source.name).read_bytes(), source.read_bytes())
            self.assertEqual((assets / manifest.name).read_bytes(), manifest.read_bytes())
            self.assertIn("Retained release assets:", result.stderr)
            recorded = [json.loads(line) for line in calls.read_text().splitlines()]
            self.assertEqual(sum(call[:2] == ["scripts/build_python_wheels.py", "check-source"] for call in recorded), 2)
            self.assertFalse(any("prepare" in call for call in recorded))
            self.assertFalse(any("build" in call for call in recorded))

    def test_release_assembly_requires_the_complete_wheel_matrix(self):
        script = (Path(__file__).resolve().parents[1] / "prepare-release.sh").read_text()
        self.assertNotIn("--allow-incomplete", script)
        self.assertNotIn("build_python_wheels.py prepare", script)

    def test_release_archives_exclude_macos_metadata_only_with_bsdtar(self):
        script = (Path(__file__).resolve().parents[1] / "prepare-release.sh").read_text()
        self.assertIn('case "$(tar --version 2>/dev/null || true)" in\n'
                      '  *bsdtar*) SMSD_RELEASE_TAR_ARGS=(--no-mac-metadata --no-xattrs) ;;\nesac\n', script)
        archives = [line for line in script.splitlines() if " tar " in line and "-czf" in line]
        self.assertEqual(len(archives), 2)
        for line in archives:
            self.assertTrue(line.startswith('COPYFILE_DISABLE=1 tar ${SMSD_RELEASE_TAR_ARGS[@]+"${SMSD_RELEASE_TAR_ARGS[@]}"} -czf '), line)


if __name__ == "__main__":
    unittest.main()
