# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""The manual Windows conda build: its limits, pins, and checks of the source and the packages."""
import hashlib
import io
import json
import os
from pathlib import Path
import importlib.util
import re
import runpy
import shutil
import subprocess
import sys
import sysconfig
import tarfile
import tempfile
import textwrap
import types
import unittest
from unittest import mock
import zipfile

try:
    from compression import zstd
    compress = zstd.compress
except ImportError:
    try:
        import zstandard
        compress = zstandard.ZstdCompressor().compress
    except ImportError:
        compress = None


ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = (ROOT / ".github/workflows/conda-windows.yml").read_text(encoding="utf-8")
RECIPE = (ROOT / "packaging/conda/meta.yaml").read_text(encoding="utf-8")
LAUNCHER_CHECK = ROOT / "packaging/conda/check_launcher.py"
SHA = "1234567890abcdef" * 2 + "12345678"
VERSION = "8.2.0"
TAGS = ("py312", "py313", "py314", "py315")
RUNNER_PATHS = {"USERPROFILE": "C:\\Users\\runneradmin", "GITHUB_WORKSPACE": "D:\\a\\SMSD\\SMSD",
                "RUNNER_TEMP": "D:\\a\\_temp", "RUNNER_TOOL_CACHE": "C:\\hostedtoolcache\\windows"}


def step(name):
    """The text of one workflow step, by name."""
    return WORKFLOW.split("      - name: " + name + "\n", 1)[1].split("\n      - ", 1)[0]


def snippet(name):
    """The inline Python program of one workflow step."""
    body = WORKFLOW.split("      - name: " + name + "\n", 1)[1].split("<<'PY'\n", 1)[1].split("          PY\n", 1)[0]
    return compile(textwrap.dedent(body), name, "exec")


BUILD = step("Build and test smsd for CPython 3.12, 3.13, 3.14 and 3.15")
INSTALL = step("Install the pinned Miniforge and conda-build under C:\\b")


class WorkflowLimitsTest(unittest.TestCase):
    def test_one_manual_windows_job_takes_only_the_release_ref(self):
        triggers = WORKFLOW.split("\non:\n", 1)[1].split("\n\n", 1)[0]
        self.assertEqual(re.findall(r"^  ([a-z_]+):", triggers, re.M), ["workflow_dispatch"])
        self.assertEqual(re.findall(r"^      ([a-z_]+):\n", triggers, re.M), ["release_tag"])
        self.assertEqual(re.findall(r"^  ([A-Za-z0-9_-]+):\n", WORKFLOW.split("\njobs:\n", 1)[1], re.M), ["win-64"])
        self.assertEqual(re.findall(r"runs-on: (\S+)", WORKFLOW), ["windows-2022"])
        self.assertLessEqual(int(re.search(r"timeout-minutes: (\d+)", WORKFLOW)[1]), 90)
        # The recipe is that of the dispatched branch, never of the release ref.
        self.assertNotIn("ref: ${{ inputs.release_tag }}", WORKFLOW)

    def test_a_second_request_for_a_release_waits_for_the_first(self):
        self.assertIn("\nconcurrency:\n  group: conda-windows-${{ inputs.release_tag }}\n"
                      "  cancel-in-progress: false\n", WORKFLOW)
        self.assertEqual(WORKFLOW.count("concurrency:"), 1)
        self.assertLess(WORKFLOW.index("\nconcurrency:"), WORKFLOW.index("\njobs:"))

    def test_the_checkout_holds_the_history_needed_to_find_master(self):
        checkout = WORKFLOW.split("uses: actions/checkout@", 1)[1].split("\n      - ", 1)[0]
        self.assertIn("\n          fetch-depth: 0\n", checkout)
        self.assertTrue(checkout.endswith("\n          persist-credentials: false"), checkout)

    def test_no_secrets_publishing_or_unreviewed_actions(self):
        lowered = WORKFLOW.lower()
        for marker in ("secrets.", "github.token", "gh_token", "id-token", "anaconda", "twine", "pypi",
                       "gh release", "gh api", "upload-file", "docker", "mvn"):
            self.assertNotIn(marker, lowered)
        self.assertNotRegex(WORKFLOW, r"(?m)^    environment:")
        self.assertEqual(set(re.findall(r"uses: ([\w.-]+/[\w.-]+)@", WORKFLOW)),
                         {"actions/checkout", "actions/setup-python", "actions/upload-artifact"})
        self.assertEqual(re.findall(r"^          path: (.+)$", WORKFLOW, re.M), ["conda-win-64/", "conda-win-64/logs/"])

    def test_builds_only_the_cpu_recipe_for_each_python_from_conda_forge(self):
        self.assertNotIn("cuda", WORKFLOW.lower())
        self.assertIn("Get-ChildItem packaging\\conda -File | Copy-Item -Destination C:\\b\\recipe", BUILD)
        self.assertIn("if (Get-ChildItem C:\\b\\recipe -Directory) { throw", BUILD)
        self.assertIn("foreach ($python in '3.12', '3.13', '3.14', '3.15')", BUILD)
        commands = re.findall(r"conda\.exe build (?!purge)(.*?) >", WORKFLOW)
        self.assertEqual(commands, ["C:\\b\\recipe -m C:\\b\\smsd-conda-variants.yaml --python $python "
                                    "--croot C:\\b\\croot --output-folder C:\\b\\channel "
                                    "-c conda-forge --override-channels"])
        for command in re.findall(r"conda\.exe (?:build (?!purge)|install ).*", WORKFLOW):
            self.assertEqual(re.findall(r" -c (\S+)", command), ["conda-forge"], command)
            self.assertIn("--override-channels", command)
        self.assertIn("Select-String -Path $log -SimpleMatch 'TEST END:'", BUILD)

    def test_installer_and_build_tools_are_pinned_and_checked(self):
        self.assertRegex(INSTALL, r"\$url = 'https://github\.com/conda-forge/miniforge/releases/download/"
                                  r"\d+\.\d+\.\d+-\d+/Miniforge3-Windows-x86_64\.exe'\n")
        self.assertRegex(INSTALL, r"\$sha = '[0-9a-f]{64}'\n")
        self.assertLess(INSTALL.index("Get-FileHash $installer"), INSTALL.index("Start-Process"))
        self.assertIn("throw 'Miniforge installer checksum differs'", INSTALL)
        self.assertRegex(INSTALL, r"conda\.exe install -y -n base -c conda-forge --override-channels "
                                  r"conda-build==\d+\.\d+\.\d+ >")

    def test_every_recorded_build_path_is_under_the_neutral_root(self):
        self.assertIn("CONDA_PKGS_DIRS: C:\\b\\pkgs\n", WORKFLOW)
        self.assertIn("'/D=C:\\b\\mf'", INSTALL)
        self.assertIn("$env:TEMP = 'C:\\b\\tmp'", BUILD)
        self.assertIn("$env:TMP = 'C:\\b\\tmp'", BUILD)
        self.assertIn("Where-Object { $_ -notmatch 'hostedtoolcache' }", BUILD)
        self.assertLess(BUILD.index("$env:PATH"), BUILD.index("conda.exe build"))
        self.assertLess(BUILD.index("$env:TMP"), BUILD.index("conda.exe build"))
        source = step("Download, check and configure the frozen source archive")
        self.assertIn('url="https://github.com/asad/SMSD/archive/$SMSD_COMMIT.tar.gz"', source)
        self.assertIn("python packaging/conda/configure_source.py", source)
        self.assertIn("--output C:/b/smsd-conda-variants.yaml", source)
        self.assertLess(source.index("<<'PY'"), source.index("configure_source.py"))


class ReleaseRefTest(unittest.TestCase):
    CODE = snippet("Validate the release ref")

    ANCESTRY = ["git", "merge-base", "--is-ancestor", SHA, "refs/remotes/origin/master"]

    def run_check(self, ref, resolved=SHA, tagged=SHA, ancestry=0):
        with tempfile.TemporaryDirectory() as directory:
            environment = Path(directory) / "github-env"
            with mock.patch.dict(os.environ, {"RELEASE_REF": ref, "GITHUB_ENV": str(environment)}), \
                    mock.patch.object(subprocess, "check_output", side_effect=[resolved + "\n", tagged + "\n"]) as git, \
                    mock.patch.object(subprocess, "run", return_value=subprocess.CompletedProcess([], ancestry)) as run:
                try:
                    exec(self.CODE, {})
                finally:
                    self.calls, self.runs = git.call_args_list, run.call_args_list
            return environment.read_text()

    def test_a_release_tag_names_its_commit_whether_or_not_master_holds_it(self):
        self.assertEqual(self.run_check("v8.2.0"), "SMSD_COMMIT=" + SHA + "\n")
        self.assertEqual([call.args[0] for call in self.calls],
                         [["git", "rev-parse", "--verify", "v8.2.0^{commit}"],
                          ["git", "rev-parse", "--verify", "refs/tags/v8.2.0^{commit}"]])
        self.assertEqual(self.runs, [])

    def test_a_full_sha_must_be_on_master(self):
        self.assertEqual(self.run_check(SHA), "SMSD_COMMIT=" + SHA + "\n")
        self.assertEqual([call.args[0] for call in self.calls], [["git", "rev-parse", "--verify", SHA + "^{commit}"]])
        self.assertEqual([call.args[0] for call in self.runs], [self.ANCESTRY])

    def test_branches_short_and_pre_release_refs_are_refused_before_git_runs(self):
        for ref in ("master", "v8.2", "v8.2.0-rc1", SHA[:12], "v8.2.0\nX"):
            with self.subTest(ref=ref), self.assertRaisesRegex(SystemExit, "stable release tag or exact"):
                self.run_check(ref)
            self.assertEqual(self.calls, [])
            self.assertEqual(self.runs, [])

    def test_a_version_name_that_is_not_a_tag_here_is_refused(self):
        with self.assertRaisesRegex(SystemExit, "not a release tag of this repository"):
            self.run_check("v8.2.0", tagged="0" * 40)
        self.assertEqual(self.runs, [])

    def test_only_commits_on_master_are_built(self):
        for ancestry, message in ((1, "The release commit is not on master"),
                                  (128, "Could not compare the release commit with master")):
            with self.subTest(ancestry=ancestry):
                with self.assertRaisesRegex(SystemExit, message):
                    self.run_check(SHA, ancestry=ancestry)
                self.assertEqual([call.args[0] for call in self.runs], [self.ANCESTRY])


class SourceArchiveTest(unittest.TestCase):
    CODE = snippet("Download, check and configure the frozen source archive")
    FILES = {"pyproject.toml": b'version = "8.2.0"\n', "packaging/conda/meta.yaml": b"package: {}\n",
             "scripts/run.bat": b"@echo off\r\n"}

    @staticmethod
    def blob(content):
        return hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest()

    def run_check(self, files=None, tree=None, comment=SHA, top="SMSD-" + SHA, link=False):
        files = self.FILES if files is None else files
        tree = self.FILES if tree is None else tree
        listing = b"".join(b"100644 blob " + self.blob(content).encode() + b"\t" + name.encode() + b"\0"
                           for name, content in tree.items())
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory) / (SHA + ".tar.gz")
            with tarfile.open(archive, "w:gz", format=tarfile.PAX_FORMAT, pax_headers={"comment": comment}) as output:
                folder = tarfile.TarInfo(top)
                folder.type = tarfile.DIRTYPE
                output.addfile(folder)
                for name, content in files.items():
                    member = tarfile.TarInfo(top + "/" + name)
                    member.size = len(content)
                    output.addfile(member, io.BytesIO(content))
                if link:
                    member = tarfile.TarInfo(top + "/link")
                    member.type, member.linkname = tarfile.SYMTYPE, "pyproject.toml"
                    output.addfile(member)
            with mock.patch.object(sys, "argv", ["-", str(archive), SHA]), \
                    mock.patch.object(subprocess, "check_output", return_value=listing) as git, \
                    mock.patch("builtins.print"):
                exec(self.CODE, {})
            git.assert_called_once_with(["git", "ls-tree", "-r", "-z", "--full-tree", SHA])

    def test_the_archive_must_hold_exactly_the_files_of_the_commit(self):
        self.run_check()
        changed = {**self.FILES, "pyproject.toml": b'version = "8.2.1"\n'}
        for case in ({"files": changed}, {"files": {**self.FILES, "extra.txt": b"x"}},
                     {"tree": {**self.FILES, "missing.txt": b"x"}},
                     {"files": {**self.FILES, "scripts/run.bat": b"@echo off\n"}}):
            with self.subTest(case=list(case)), self.assertRaisesRegex(SystemExit, "Archive differs from the release commit"):
                self.run_check(**case)

    def test_the_archive_must_name_the_commit_and_hold_only_files(self):
        with self.assertRaisesRegex(SystemExit, "names another commit"):
            self.run_check(comment="0" * 40)
        with self.assertRaisesRegex(SystemExit, "Unexpected archive path"):
            self.run_check(top="SMSD-master")
        with self.assertRaisesRegex(SystemExit, "Unexpected archive member"):
            self.run_check(link=True)


@unittest.skipIf(compress is None, "needs a Zstandard compressor")
class PackageScanTest(unittest.TestCase):
    CODE = snippet("Check, scan and checksum the packages")

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.channel = self.root / "channel"
        self.work = self.root / "work"
        self.work.mkdir()
        (self.root / "smsd-conda-variants.yaml").write_text(json.dumps(
            {"python": ["3.12", "3.13", "3.14", "3.15"], "smsd_version": [VERSION], "smsd_source_commit": [SHA]}))
        self.reset()

    def tearDown(self):
        self.directory.cleanup()

    def reset(self):
        shutil.rmtree(self.channel, ignore_errors=True)
        for tag in TAGS:
            self.package(tag)

    @staticmethod
    def tar(files):
        stream = io.BytesIO()
        with tarfile.open(fileobj=stream, mode="w") as archive:
            for name, content in files.items():
                member = tarfile.TarInfo(name)
                member.size = len(content)
                archive.addfile(member, io.BytesIO(content))
        return compress(stream.getvalue())

    def package(self, tag, payload=None, subdir="win-64", name=None, index=None, commit=SHA, paths=(),
                info=None, paths_json=True):
        minor = int(tag[3:])
        stem = name or f"smsd-{VERSION}-{tag}h1a2b3c4_0"
        metadata = {"name": "smsd", "version": VERSION, "subdir": "win-64", "build_number": 0,
                    "depends": [f"python >=3.{minor},<3.{minor + 1}.0a0", f"python_abi 3.{minor}.* *_cp3{minor}"]}
        info = {"info/index.json": json.dumps({**metadata, **(index or {})}).encode(),
                "info/about.json": json.dumps({"extra": {"source_commit": commit}}).encode(),
                **({"info/paths.json": json.dumps({"paths": list(paths)}).encode()} if paths_json else {}),
                **(info or {})}
        files = {"Lib/site-packages/smsd/__init__.py": b"__version__ = '8.2.0'\n", **(payload or {})}
        path = self.channel / subdir / (stem + ".conda")
        path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path, "w") as output:
            output.writestr("metadata.json", '{"conda_pkg_format_version": 2}')
            output.writestr(f"info-{stem}.tar.zst", self.tar(info))
            output.writestr(f"pkg-{stem}.tar.zst", self.tar(files))
        return path

    def run_check(self, runner=RUNNER_PATHS):
        original = Path.cwd()
        try:
            os.chdir(self.work)
            with mock.patch.dict(os.environ, {"SMSD_CONDA_ROOT": str(self.root), "SMSD_COMMIT": SHA}), \
                    mock.patch("builtins.print"):
                for variable in RUNNER_PATHS:
                    os.environ.pop(variable, None)
                os.environ.update(runner)
                exec(self.CODE, {})
        finally:
            os.chdir(original)

    def assertRefused(self, message, **options):
        with self.assertRaisesRegex(SystemExit, message):
            self.run_check(**options)
        self.assertFalse((self.work / "conda-win-64").exists())

    def test_clean_packages_are_copied_with_their_checksums(self):
        self.run_check()
        output = self.work / "conda-win-64"
        sums = (output / "SHA256SUMS").read_text().splitlines()
        self.assertEqual(len(sums), 4)
        for line, tag in zip(sums, TAGS):
            digest, name = line.split("  ")
            self.assertEqual(name, f"win-64/smsd-{VERSION}-{tag}h1a2b3c4_0.conda")
            self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), digest)
            self.assertEqual((output / name).read_bytes(), (self.channel / name).read_bytes())
        self.assertEqual((output / "smsd-conda-variants.yaml").read_text(),
                         (self.root / "smsd-conda-variants.yaml").read_text())

    def test_runner_user_and_workspace_paths_are_refused(self):
        payloads = ({"Lib/site-packages/smsd/a.pyc": b"\0C:\\Users\\runneradmin\\AppData\\Local\\Temp\\x.py"},
                    {"Lib/site-packages/smsd/_smsd.pyd": "D:\\a\\SMSD\\SMSD\\cpp\\x.hpp".encode("utf-16-le")},
                    {"info/x.json": b'{"cwd": "c:\\\\users\\\\someone\\\\src"}'},
                    {"Lib/site-packages/smsd/b.txt": b"see c:/Users/x"},
                    {"Lib/site-packages/smsd/c.txt": b"C:\\hostedtoolcache\\windows\\Python"},
                    {"Lib/site-packages/runneradmin/d.txt": b""},
                    {"Lib/site-packages/smsd/e.sh": b"cd /c/Users/someone/src"},
                    {"Lib/site-packages/smsd/f.sh": "cd /d/a/SMSD/SMSD/cpp".encode("utf-16-le")})
        for payload in payloads:
            with self.subTest(payload=list(payload)):
                self.package("py313", payload)
                self.assertRefused("Package check failed")
        self.package("py313")
        self.run_check()

    def test_this_runners_own_paths_are_refused_wherever_they_lie(self):
        runner = {"USERPROFILE": "F:\\Home\\builder", "GITHUB_WORKSPACE": "/e/w/SMSD",
                  "RUNNER_TEMP": "E:\\w\\_temp\\", "RUNNER_TOOL_CACHE": "G:/tools"}
        for path in ("F:\\Home\\builder\\x", "e:\\w\\smsd\\cpp\\x.hpp", "E:/w/SMSD/x", "E:\\\\w\\\\SMSD\\\\x",
                     "E:/w/_temp/y", "G:\\tools\\Python"):
            for encoding in ("utf-8", "utf-16-le"):
                with self.subTest(path=path, encoding=encoding):
                    self.package("py313", {"Lib/site-packages/smsd/a.txt": path.encode(encoding)})
                    self.assertRefused("Package check failed", runner=runner)
        self.package("py313", {"Lib/site-packages/smsd/a.txt": b"E:\\w\\other F:\\Home\\someone G:\\Python"})
        self.run_check(runner=runner)
        for value in ("C:\\", "/c/", "_temp", "/home/runner/work"):
            with self.subTest(value=value), self.assertRaisesRegex(SystemExit, "RUNNER_TEMP does not name a folder"):
                self.run_check(runner={**RUNNER_PATHS, "RUNNER_TEMP": value})

    BUILD = "C:\\b\\croot\\smsd_1793000000000\\_h_env"
    LAUNCHER = b"MZ\0launcher#!" + (BUILD + "\\python.exe").encode() + b"\r\n"

    def test_build_paths_outside_relocated_files_are_refused(self):
        payloads = ({"Scripts/smsd-mcp.exe": self.LAUNCHER},
                    {"Lib/site-packages/smsd/_smsd.pyd": "C:\\b\\croot\\smsd_1\\work\\cpp\\x.hpp".encode("utf-16-le")},
                    {"Lib/site-packages/smsd/a.json": b'{"root": "c:\\\\b\\\\croot"}'},
                    {"Lib/site-packages/smsd/b.txt": b"C:/b/recipe/build.py"},
                    {"Lib/site-packages/other-1.0.dist-info/direct_url.json": b'{"url": "file:///C:/b/croot/x/work"}'})
        for payload in payloads:
            with self.subTest(payload=list(payload)):
                self.package("py314", payload)
                self.assertRefused(re.escape(list(payload)[0]) + " holds the build path .* outside a relocated file")
        # conda copies a listed file unchanged unless it names a placeholder.
        self.package("py314", {"Scripts/smsd-mcp.exe": self.LAUNCHER},
                     paths=[{"_path": "Scripts/smsd-mcp.exe", "path_type": "hardlink"}])
        self.assertRefused("Scripts/smsd-mcp.exe holds the build path")
        self.package("py314", paths_json=False)
        self.assertRefused("info/paths.json is missing")

    def test_the_pip_record_is_accepted_in_either_windows_layout(self):
        record = b'{"dir_info": {}, "url": "file:///C:/b/croot/smsd_1793000000000/work"}'
        for prefix in ("Lib", "lib/python"):
            with self.subTest(prefix=prefix):
                shutil.rmtree(self.work / "conda-win-64", ignore_errors=True)
                self.package("py315", {f"{prefix}/site-packages/smsd-{VERSION}.dist-info/direct_url.json": record})
                self.run_check()

    def test_relocated_and_inert_build_paths_are_accepted(self):
        script = b"#!" + (self.BUILD + "\\python.exe").encode() + b"\nfrom smsd.mcp_server import main\n"
        self.package("py314", {"Scripts/smsd-mcp.exe": self.LAUNCHER, "Scripts/smsd-mcp-script.py": script,
                               "Lib/site-packages/smsd/__pycache__/x.cpython-314.pyc": b"\0C:\\b\\croot\\x.py",
                               f"Lib/site-packages/smsd-{VERSION}.dist-info/direct_url.json":
                                   b'{"dir_info": {}, "url": "file:///C:/b/croot/smsd_1793000000000/work"}'},
                     paths=[{"_path": "Scripts/smsd-mcp.exe", "file_mode": "binary", "prefix_placeholder": self.BUILD},
                            {"_path": "Scripts/smsd-mcp-script.py", "file_mode": "text",
                             "prefix_placeholder": self.BUILD}],
                     info={"info/has_prefix": (self.BUILD + " binary Scripts/smsd-mcp.exe\n").encode(),
                           "info/recipe/meta.yaml": b"build:\n  script: C:\\b\\recipe/build.py\n"})
        self.run_check()
        self.assertTrue((self.work / "conda-win-64" / "SHA256SUMS").is_file())

    def test_only_one_cpu_package_for_each_python_is_accepted(self):
        cases = (lambda: self.package("py312", name=f"smsd-cuda-{VERSION}-py312h1a2b3c4_0"),
                 lambda: self.package("py312", subdir="noarch"),
                 lambda: self.package("py312", name=f"smsd-{VERSION}-py312h9a9b9c9_0"),
                 lambda: (self.channel / "win-64" / f"smsd-{VERSION}-py312h1a2b3c4_0.tar.bz2").write_bytes(b""))
        for case in cases:
            with self.subTest(case=cases.index(case)):
                self.reset()
                case()
                self.assertRefused("Unexpected package")
        self.reset()
        (self.channel / "win-64" / f"smsd-{VERSION}-py315h1a2b3c4_0.conda").unlink()
        self.assertRefused("Expected one smsd package for each of CPython 3.12-3.15")

    def test_package_identity_python_and_source_commit_are_checked(self):
        for case, message in (({"index": {"subdir": "linux-64"}}, "unexpected package identity"),
                              ({"index": {"build_number": 1}}, "unexpected package identity"),
                              ({"index": {"depends": ["python >=3.12,<3.13.0a0"]}}, "expected the dependency"),
                              ({"commit": "0" * 40}, "recorded source commit differs")):
            with self.subTest(case=case):
                self.package("py314", **case)
                self.assertRefused(message)
        self.package("py314")
        variants = self.root / "smsd-conda-variants.yaml"
        variants.write_text(variants.read_text().replace(SHA, "0" * 40))
        self.assertRefused("variant inputs name another source commit")


class RecipeLauncherTest(unittest.TestCase):
    MODULE = "/env/Lib/site-packages/smsd/mcp_server.py"
    MISSING_SDK = "ImportError: The MCP server requires the official SDK; install smsd[mcp]\n"

    def test_only_the_windows_recipe_test_runs_the_launcher(self):
        test = RECIPE.split("\ntest:\n", 1)[1].split("\n\n", 1)[0]
        self.assertIn("\n  files:  # [win]\n    - check_launcher.py  # [win]\n", test)
        commands = test.split("\n  commands:\n", 1)[1]
        self.assertEqual(re.findall(r"^    - (.+)$", commands, re.M),
                         ["python -m pip check", "python check_launcher.py  # [win]"])
        # Each added line carries the Windows selector, so Linux and macOS render the recipe as before.
        added = [line for line in RECIPE.splitlines() if "check_launcher" in line or line.startswith("  files:")]
        self.assertEqual(len(added), 3)
        for line in added:
            self.assertTrue(line.endswith("  # [win]"), line)
        # The workflow copies the files of packaging/conda, so the check reaches the Windows recipe.
        self.assertTrue(LAUNCHER_CHECK.is_file())

    def run_launcher_check(self, returncode=0, stderr="", present=True):
        with tempfile.TemporaryDirectory() as scripts:
            launcher = Path(scripts, "smsd-mcp.exe" if os.name == "nt" else "smsd-mcp")
            if present:
                launcher.write_bytes(b"")
            package, server = types.ModuleType("smsd"), types.ModuleType("smsd.mcp_server")
            server.__file__, package.mcp_server = self.MODULE, server
            result = subprocess.CompletedProcess([], returncode, "", stderr)
            with mock.patch.dict(sys.modules, {"smsd": package, "smsd.mcp_server": server}), \
                    mock.patch.object(sysconfig, "get_path", return_value=scripts), \
                    mock.patch.object(subprocess, "run", return_value=result) as run, mock.patch("builtins.print"):
                try:
                    runpy.run_path(str(LAUNCHER_CHECK), run_name="__main__")
                finally:
                    self.runs = [call.args[0] for call in run.call_args_list]
            self.assertEqual(self.runs, [[str(launcher), "--help"]] if present else [])

    def traceback(self, module):
        return ('Traceback (most recent call last):\n  File "' + module + '", line 461, in create_server\n'
                + self.MISSING_SDK)

    def test_a_launcher_running_this_environment_passes(self):
        self.run_launcher_check(1, self.traceback(self.MODULE))

    def test_a_launcher_running_another_environment_or_failing_otherwise_is_refused(self):
        with self.assertRaisesRegex(SystemExit, "loaded SMSD from another environment"):
            self.run_launcher_check(1, self.traceback("C:\\b\\croot\\smsd_1\\_h_env\\Lib\\site-packages\\smsd\\mcp_server.py"))
        for returncode, stderr in ((0, ""), (0, self.traceback(self.MODULE)), (1, ""),
                                   (1, "Fatal error in launcher: Unable to create process using C:\\b\\croot\\python.exe\n"),
                                   (1, self.traceback(self.MODULE) + "SystemExit: 2\n")):
            with self.subTest(returncode=returncode, stderr=stderr), \
                    self.assertRaisesRegex(SystemExit, "did not stop at the missing MCP SDK"):
                self.run_launcher_check(returncode, stderr)
        with self.assertRaisesRegex(SystemExit, "is missing"):
            self.run_launcher_check(present=False)

    @unittest.skipUnless(importlib.util.find_spec("smsd") and not importlib.util.find_spec("mcp")
                         and Path(sysconfig.get_path("scripts"), "smsd-mcp.exe" if os.name == "nt" else "smsd-mcp").is_file(),
                         "needs an installed smsd-mcp launcher without the MCP SDK")
    def test_the_check_passes_against_an_installed_launcher(self):
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run([sys.executable, "-I", str(LAUNCHER_CHECK)], cwd=directory,
                                    capture_output=True, text=True, timeout=300)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("smsd-mcp", result.stdout)


if __name__ == "__main__":
    unittest.main()
