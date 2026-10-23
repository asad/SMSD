# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Release preparation must preserve an existing collection of assets."""
import os
from pathlib import Path
import shutil
import subprocess
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


if __name__ == "__main__":
    unittest.main()
