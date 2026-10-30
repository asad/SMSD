# SPDX-License-Identifier: Apache-2.0
"""CUDA collection derives success from actual case identities and installed wheel bytes."""

import importlib.util
import json
import os
from pathlib import Path
import sys
import types
import xml.etree.ElementTree as ET
import zipfile

import pytest

SPEC = importlib.util.spec_from_file_location("cuda_proof", Path(__file__).parents[1] / "cuda_wheel_provenance.py")
PROOF = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PROOF)


def write_xml(path, cases):
    root = ET.Element("testsuite", name="results", tests=str(len(cases)), failures="0", errors="0",
                      skipped=str(sum(entry[2] for entry in cases)))
    for entry in cases:
        classname, name, skipped = entry[:3]
        case = ET.SubElement(root, "testcase", classname=classname, name=name)
        if skipped:
            ET.SubElement(case, "skipped", message=entry[4] if len(entry) > 4 else "optional dependency unavailable")
        if len(entry) >= 4:
            properties = ET.SubElement(case, "properties")
            ET.SubElement(properties, "property", name="smsd_nodeid", value=entry[3])
    ET.ElementTree(root).write(path)


@pytest.fixture
def checked_results(tmp_path):
    prefix = tmp_path / "env"
    module = prefix / "smsd_cuda/__init__.py"
    module.parent.mkdir(parents=True)
    module.write_bytes(b"version fixture\n")
    library = module.with_name("libsmsd_cuda_screen.so")
    library.write_bytes(b"library fixture")
    native = prefix / "smsd/_smsd.so"
    native.parent.mkdir()
    native.write_bytes(b"core native fixture")
    core_module = native.with_name("__init__.py")
    core_module.write_bytes(b"core python fixture")
    cpu_wheel = tmp_path / "smsd-8.1.0-cp314-cp314-linux_x86_64.whl"
    with zipfile.ZipFile(cpu_wheel, "w") as archive:
        archive.writestr("smsd/__init__.py", core_module.read_bytes())
        archive.writestr("smsd/_smsd.so", native.read_bytes())
    wheel = tmp_path / "smsd_cuda-8.1.0-py3-none-manylinux_2_34_x86_64.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("smsd_cuda/__init__.py", module.read_bytes())
        archive.writestr("smsd_cuda/libsmsd_cuda_screen.so", library.read_bytes())
    source = {"source_commit": "a" * 40, "sdist_sha256": "b" * 64, "version": "8.1.0"}
    sdk, checkout = tmp_path / "sdk", tmp_path / "checkout"
    core_tests, addon_tests = sdk / "python/tests", checkout / "packaging/smsd-cuda/tests"
    core_tests.mkdir(parents=True)
    addon_tests.mkdir(parents=True)
    modules = [core_tests / "test_smsd.py", core_tests / "test_optional.py",
               addon_tests / "test_smsd_cuda.py", addon_tests / "test_resident_fingerprints.py"]
    for path in modules:
        path.write_text("fixture test source\n")
    relative = [os.path.relpath(path, sdk) for path in modules]
    nodeids = [relative[0] + "::TestParsing::test_parse", relative[0] + "::TestParsing::test_other",
               relative[2] + "::test_smsd_attaches_the_device",
               relative[3] + "::test_resident_integer_thresholds_and_order[32-5]"]
    native_xml, python_xml = tmp_path / "native.xml", tmp_path / "python.xml"
    write_xml(native_xml, [(name, name, False) for name in sorted(PROOF.NATIVE_TESTS)])
    write_xml(python_xml, [("python.tests.test_smsd.TestParsing", "test_parse", False, nodeids[0]),
                           ("python.tests.test_smsd.TestParsing", "test_other", False, nodeids[1]),
                           ("tests.test_smsd_cuda", "test_smsd_attaches_the_device", False, nodeids[2]),
                           ("tests.test_resident_fingerprints", "test_resident_integer_thresholds_and_order[32-5]", False, nodeids[3]),
                           ("", "python.tests.test_optional", True)])
    native_dir = tmp_path / "native-build"
    native_dir.mkdir()
    for name in PROOF.NATIVE_TESTS:
        (native_dir / name).write_bytes(b"native executable fixture")
    note = {"source_commit": source["source_commit"], "sdist_sha256": source["sdist_sha256"],
            "cpu_wheel": str(cpu_wheel), "cpu_wheel_sha256": PROOF.sha256(cpu_wheel),
            "cuda_wheel_sha256": PROOF.sha256(wheel), "sdk_source": str(sdk), "checkout": str(checkout),
            "native_dir": str(native_dir), "native_binaries": {name: PROOF.sha256(native_dir / name) for name in PROOF.NATIVE_TESTS},
            "cuda_toolkit": "Cuda compilation tools, release 12.8, V12.8.93", "finished_ns": 0}
    note_path = tmp_path / "build-inputs.json"
    note_path.write_text(json.dumps(note))
    runtime_path = prefix / "libcudart.so.12"
    runtime_path.write_bytes(b"runtime fixture")
    identity = {"source_commit": source["source_commit"], "sdist_sha256": source["sdist_sha256"],
                "smsd_version": source["version"], "wheel": wheel.name, "wheel_sha256": PROOF.sha256(wheel),
                "prefix": str(prefix), "module": str(module), "module_sha256": PROOF.sha256(module),
                "library": str(library), "library_sha256": PROOF.sha256(library),
                "native_module": str(native), "native_module_sha256": PROOF.sha256(native),
                "pytest_exit_status": 0, "junit_xml": str(python_xml), "core_module": str(core_module),
                "cpu_wheel": str(cpu_wheel), "cpu_wheel_sha256": PROOF.sha256(cpu_wheel),
                "core_files": PROOF.installed_files(cpu_wheel, "smsd", prefix, core_module, native),
                "cuda_files": PROOF.installed_files(wheel, "smsd_cuda", prefix, module, library),
                "build_note": str(note_path), "build_note_sha256": PROOF.sha256(note_path),
                "native_xml": str(native_xml), "native_xml_sha256": PROOF.sha256(native_xml),
                "ctest_exit_status": 0, "ctest_started_ns": 1, "ctest_finished_ns": native_xml.stat().st_mtime_ns + 1,
                "runtime": {"version": 12080, "package_version": "12.8.90", "libraries": {str(runtime_path): PROOF.sha256(runtime_path)}},
                "test_roots": [str(core_tests), str(addon_tests)],
                "collection": {"items": nodeids.copy(), "selected": nodeids.copy(), "deselected": [],
                               "modules": {name: {"outcome": "skipped" if index == 1 else "passed", "reason": "optional dependency unavailable" if index == 1 else ""}
                                           for index, name in enumerate(relative)},
                               "outcomes": {nodeid: {"outcome": "passed", "reason": ""} for nodeid in nodeids},
                               "rootdir": str(sdk), "filters": {}, "env_addopts": "", "ini_addopts": ["--import-mode=importlib"],
                               "selectors": {"python_files": ["test_*.py", "*_test.py"], "python_classes": ["Test"], "python_functions": ["test_*"]},
                               "args": [str(core_tests), str(addon_tests)],
                               "invocation_args": [str(core_tests), str(addon_tests), "-q", "-rs", "-p", "no:cacheprovider", "-p", "cuda_wheel_provenance",
                                                   "--junitxml=" + str(python_xml), "--rootdir", str(sdk)]}}
    config = sdk / "pyproject.toml"
    config.write_text('[tool.pytest.ini_options]\naddopts = "--import-mode=importlib"\ntimeout = 60\n')
    collection = identity["collection"]
    collection["selectors"]["python_functions"] = ["test"]
    collection["invocation_args"][-2:-2] = ["-c", str(config)]
    collection["module_paths"] = {alias: str((sdk / alias).resolve()) for alias in relative}
    collection["config"] = {"config_path": str(config), "config_sha256": PROOF.sha256(config),
        "importmode": "importlib", "timeout": "60", "ini_addopts": ["--import-mode=importlib"],
        "invocation_args": collection["invocation_args"].copy(), "rootdir": str(sdk)}
    collection["pytest_version"] = "9.1.1"
    collection["pytest_config_source"] = {"path": str(config), "sha256": PROOF.sha256(config)}
    identity["capabilities"] = {"available": True, "device_info": "Fixture GPU (compute 7.5, 1024 MB)",
        "binary_mma_available": 0, "binary_mma_error": "", "testing_failure_hook_present": False}
    return source, wheel, identity, native_xml, python_xml


def test_producer_records_skips_and_only_public_proof_fields(checked_results):
    proof, summary = PROOF.proof_from_results(*checked_results)
    assert proof["installed_wheel_tests"] == proof["native_gpu_tests"] == "passed"
    assert summary["native_passed"] == 6
    assert summary["python_passed"] == 4
    assert summary["python_skipped"] == [{"classname": "", "name": "python.tests.test_optional",
                                          "reason": "optional dependency unavailable"}]
    assert set(proof) == {"schema_version", "source_commit", "sdist_sha256", "smsd_version",
                          "wheel", "wheel_sha256", "installed_wheel_tests", "native_gpu_tests"}


@pytest.mark.parametrize("field, value", [("source_commit", "c" * 40), ("sdist_sha256", "c" * 64),
                                           ("wheel_sha256", "c" * 64), ("smsd_version", "8.0.3"),
                                           ("pytest_exit_status", 1), ("pytest_exit_status", False)])
def test_reject_changed_invocation_identity(checked_results, field, value):
    checked_results[2][field] = value
    with pytest.raises(ValueError, match="identity"):
        PROOF.proof_from_results(*checked_results)


def test_reject_foreign_loaded_origin_and_changed_installed_library(checked_results, tmp_path):
    identity = checked_results[2]
    outside = tmp_path / "foreign.so"
    outside.write_bytes(Path(identity["library"]).read_bytes())
    original = identity["library"]
    identity["library"] = str(outside)
    with pytest.raises(ValueError, match="installed environment"):
        PROOF.proof_from_results(*checked_results)
    identity["library"] = original
    Path(original).write_bytes(b"changed library")
    with pytest.raises(ValueError, match="differs"):
        PROOF.proof_from_results(*checked_results)


@pytest.mark.parametrize("broken", ["failure", "error", "empty", "missing_native", "missing_device",
                                    "missing_resident", "unrecorded_skip"])
def test_reject_missing_or_failed_actual_cases(checked_results, broken):
    native_xml, python_xml = checked_results[3:]
    path = python_xml if broken in ("missing_device", "missing_resident", "unrecorded_skip") else native_xml
    tree = ET.parse(path)
    root = tree.getroot()
    if broken in ("failure", "error"):
        ET.SubElement(root.find("testcase"), broken)
    elif broken == "empty":
        root.clear()
        root.set("name", "empty")
    elif broken == "unrecorded_skip":
        root.set("skipped", "2")
    else:
        cases = root.findall("testcase")
        remove = cases[0] if broken == "missing_native" else cases[2 if broken == "missing_device" else 3]
        root.remove(remove)
        root.set("tests", str(len(cases) - 1))
    tree.write(path)
    with pytest.raises(ValueError):
        PROOF.proof_from_results(*checked_results)


def test_session_hook_captures_loaded_origin_in_test_process(checked_results, tmp_path, monkeypatch):
    source, wheel, identity, _, python_xml = checked_results
    archive = tmp_path / "smsd-8.1.0.tar.gz"
    archive.write_bytes(b"source archive fixture")
    manifest = tmp_path / "SOURCE_MANIFEST.json"
    manifest.write_text(json.dumps({**source, "sdist": archive.name}))
    output = tmp_path / "identity.json"
    monkeypatch.setenv("SMSD_CUDA_SOURCE_MANIFEST", str(manifest))
    monkeypatch.setenv("SMSD_CUDA_TEST_WHEEL", str(wheel))
    monkeypatch.setenv("SMSD_CUDA_TEST_IDENTITY", str(output))
    monkeypatch.setenv("SMSD_CUDA_BUILD_NOTE", identity["build_note"])
    monkeypatch.setenv("SMSD_CUDA_NATIVE_XML", str(checked_results[3]))
    monkeypatch.setenv("SMSD_CUDA_CTEST_STATUS", "0")
    monkeypatch.setenv("SMSD_CUDA_CTEST_STARTED_NS", "1")
    monkeypatch.setenv("SMSD_CUDA_CTEST_FINISHED_NS", str(identity["ctest_finished_ns"]))
    monkeypatch.setattr(PROOF, "_COLLECTION", {**identity["collection"], "items": set(identity["collection"]["items"])})
    monkeypatch.setattr(PROOF, "loaded_runtime", lambda cuda: identity["runtime"])
    monkeypatch.setattr(PROOF, "loaded_capabilities", lambda cuda: identity["capabilities"])
    native = types.SimpleNamespace(__version__="8.1.0", __file__=identity["native_module"])
    monkeypatch.setitem(sys.modules, "smsd", types.SimpleNamespace(__version__="8.1.0", __file__=identity["core_module"], _smsd=native))
    monkeypatch.setitem(sys.modules, "smsd_cuda", types.SimpleNamespace(
        __version__="8.1.0", __smsd_version__="8.1.0", __file__=identity["module"],
        _library=types.SimpleNamespace(_name=identity["library"])))
    monkeypatch.setattr(sys, "prefix", identity["prefix"])
    ini = {"addopts": ["--import-mode=importlib"], "timeout": 60, **identity["collection"]["selectors"]}
    session = types.SimpleNamespace(config=types.SimpleNamespace(
        inipath=Path(identity["collection"]["config"]["config_path"]),
        option=types.SimpleNamespace(xmlpath=str(python_xml), importmode="importlib"), rootpath=Path(identity["collection"]["rootdir"]),
        invocation_params=types.SimpleNamespace(args=identity["collection"]["invocation_args"]),
        args=identity["collection"]["args"], getini=lambda name: ini[name]))
    PROOF.pytest_sessionfinish(session, 0)
    recorded = PROOF.read_json(output)
    assert recorded["library"] == identity["library"]
    assert recorded["wheel_sha256"] == PROOF.sha256(wheel)
    assert recorded["sdist_sha256"] == PROOF.sha256(archive)


@pytest.mark.parametrize("change", ["deselected", "keyword", "truncated", "xml_node", "required_skip", "module_skip"])
def test_reject_incomplete_collection_or_xml_join(checked_results, change):
    identity = checked_results[2]
    collection = identity["collection"]
    if change == "deselected":
        collection["deselected"].append(collection["selected"][0])
    elif change == "keyword":
        collection["filters"]["keyword"] = "selected_case"
    elif change == "truncated":
        collection["outcomes"].pop(collection["selected"][0])
    elif change == "module_skip":
        nodeid = next(name for name in collection["modules"] if name.endswith("test_smsd_cuda.py"))
        collection["modules"][nodeid] = {"outcome": "skipped", "reason": "device unavailable"}
    else:
        tree = ET.parse(checked_results[4])
        case = tree.getroot().findall("testcase")[2]
        if change == "xml_node":
            case.find("./properties/property").set("value", "uncollected.py::test_other")
        else:
            ET.SubElement(case, "skipped", message="device unavailable")
            tree.getroot().set("skipped", "2")
            collection["outcomes"][collection["selected"][2]] = {"outcome": "skipped", "reason": "device unavailable"}
        tree.write(checked_results[4])
    with pytest.raises(ValueError):
        PROOF.proof_from_results(*checked_results)


@pytest.mark.parametrize("change", ["core_python", "core_native", "cpu_wheel", "native_binary", "stale_ctest", "runtime_floor"])
def test_reject_changed_build_or_installed_inputs(checked_results, change):
    identity = checked_results[2]
    if change in ("core_python", "core_native", "cpu_wheel"):
        key = {"core_python": "core_module", "core_native": "native_module", "cpu_wheel": "cpu_wheel"}[change]
        Path(identity[key]).write_bytes(b"changed input")
    elif change == "native_binary":
        note = PROOF.read_json(identity["build_note"])
        (Path(note["native_dir"]) / sorted(PROOF.NATIVE_TESTS)[0]).write_bytes(b"changed binary")
    elif change == "stale_ctest":
        identity["ctest_started_ns"] = checked_results[3].stat().st_mtime_ns + 1
    else:
        identity["runtime"]["version"] = 12090
    with pytest.raises(ValueError):
        PROOF.proof_from_results(*checked_results)


def test_reject_private_identity_inside_release_directory(checked_results, tmp_path, monkeypatch):
    release = tmp_path / "release"
    release.mkdir()
    identity = release / "installed-identity.json"
    identity.write_text(json.dumps(checked_results[2]))
    monkeypatch.setattr(sys, "argv", ["cuda_wheel_provenance.py", "--release-dir", str(release),
                                    "--wheel", str(checked_results[1]), "--identity", str(identity),
                                    "--native-xml", str(checked_results[3]), "--python-xml", str(checked_results[4]),
                                    "--evidence-summary", str(tmp_path / "summary.json")])
    with pytest.raises(ValueError, match="outside the release"):
        PROOF.main()
    assert {path.name for path in release.iterdir()} == {identity.name}


MMA_REASON = "binary MMA prototype needs an enabled build and an SM80-SM90 device"
INJECTION_REASON = "failure injection requires the separate SMSD_CUDA_TESTING library"
CAPABILITY_CASES = [("test_opt_in_binary_mma_integer_intersections[" + str(words) + "]", MMA_REASON)
                    for words in (1, 3, 4, 5, 33)] + [
    (function + "[" + str(failure) + "]", INJECTION_REASON)
    for function in ("test_resident_top_k_injected_failure_keeps_payload_and_handle",
                     "test_resident_graph_failure_recovery_preserves_payload") for failure in (1, 2)]


def add_cuda_skip(checked_results, name, reason):
    identity, xml = checked_results[2], checked_results[4]
    collection = identity["collection"]
    alias = next(alias for alias in collection["modules"] if alias.endswith("test_resident_fingerprints.py"))
    nodeid = alias + "::" + name
    for key in ("items", "selected"):
        collection[key].append(nodeid)
    collection["outcomes"][nodeid] = {"outcome": "skipped", "reason": repr((str(Path(identity["test_roots"][1]) / "test_resident_fingerprints.py"), 305, "Skipped: " + reason))}
    tree = ET.parse(xml)
    root = tree.getroot()
    root.set("tests", str(int(root.get("tests")) + 1))
    root.set("skipped", str(int(root.get("skipped")) + 1))
    case = ET.SubElement(root, "testcase", classname="test_resident_fingerprints", name=name)
    ET.SubElement(case, "skipped", message=reason)
    ET.SubElement(ET.SubElement(case, "properties"), "property", name="smsd_nodeid", value=nodeid)
    tree.write(xml)


@pytest.mark.parametrize("name,reason", CAPABILITY_CASES)
def test_only_exact_unsupported_cuda_controls_are_retained_as_skips(checked_results, name, reason):
    add_cuda_skip(checked_results, name, reason)
    proof, summary = PROOF.proof_from_results(*checked_results)
    assert proof["installed_wheel_tests"] == "passed"
    assert any(case["name"] == name and case["reason"] == reason for case in summary["python_skipped"])


@pytest.mark.parametrize("name,reason,change", [
    ("test_opt_in_binary_mma_integer_intersections[2]", MMA_REASON, None),
    ("test_opt_in_binary_mma_integer_intersections[1]", "device unavailable", None),
    ("test_opt_in_binary_mma_integer_intersections[1]", MMA_REASON, "mma_available"),
    ("test_resident_top_k_injected_failure_keeps_payload_and_handle[1]", INJECTION_REASON, "hook_present"),
    ("test_resident_graph_failure_recovery_preserves_payload[3]", INJECTION_REASON, None),
    ("test_resident_integer_thresholds_and_order[32-5]", MMA_REASON, None),
])
def test_supported_required_or_mislabelled_cuda_skips_fail(checked_results, name, reason, change):
    add_cuda_skip(checked_results, name, reason)
    if change == "mma_available":
        checked_results[2]["capabilities"]["binary_mma_available"] = 1
    elif change == "hook_present":
        checked_results[2]["capabilities"]["testing_failure_hook_present"] = True
    with pytest.raises(ValueError):
        PROOF.proof_from_results(*checked_results)


def test_portable_fallback_skip_requires_mma_capability(checked_results):
    name = "test_binary_mma_request_keeps_unsupported_device_on_portable_path"
    add_cuda_skip(checked_results, name, "this control requires the portable fallback device or build")
    with pytest.raises(ValueError):
        PROOF.proof_from_results(*checked_results)
    checked_results[2]["capabilities"]["binary_mma_available"] = 1
    assert PROOF.proof_from_results(*checked_results)[0]["installed_wheel_tests"] == "passed"


def test_actual_pytest_configuration_keeps_external_modules_and_default_test_prefix(tmp_path, monkeypatch):
    sdk, external = tmp_path / "sdk", tmp_path / "external"
    (sdk / "python/tests").mkdir(parents=True)
    external.mkdir()
    config = sdk / "pyproject.toml"
    config.write_text('[tool.pytest.ini_options]\naddopts="--import-mode=importlib"\ntimeout=60\n')
    (sdk / "python/tests/test_core.py").write_text("def testordinary(): pass\n")
    (external / "test_external.py").write_text("def test_other(): pass\n")
    captured = {}
    class Capture:
        def pytest_collectstart(self, collector):
            PROOF.pytest_collectstart(collector)

        def pytest_collection_finish(self, session):
            captured["config"] = session.config.inipath
            captured["selectors"] = session.config.getini("python_functions")
            captured["items"] = [(item.nodeid, str(item.path.resolve())) for item in session.items]
    monkeypatch.setenv("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")
    monkeypatch.setenv("SMSD_CUDA_TEST_IDENTITY", "recording-fixture")
    monkeypatch.setattr(PROOF, "_COLLECTION", {"module_paths": {}})
    status = pytest.main([str(sdk / "python/tests"), str(external), "-c", str(config),
                          "--rootdir", str(sdk), "--collect-only", "-q", "-p", "no:cacheprovider"], plugins=[Capture()])
    assert status == 0 and captured["config"] == config
    assert captured["selectors"] == ["test"]
    assert {path for _, path in captured["items"]} == {str((sdk / "python/tests/test_core.py").resolve()), str((external / "test_external.py").resolve())}
    assert any(node.endswith("::testordinary") for node, _ in captured["items"])
    assert set(PROOF._COLLECTION["module_paths"].values()) == {path for _, path in captured["items"]}


def test_actual_bare_external_nodeids_use_recorded_paths(checked_results):
    identity, xml = checked_results[2], checked_results[4]
    collection = identity["collection"]
    mapping = {alias: Path(alias).name if alias.startswith("..") else alias for alias in collection["modules"]}
    def rename(node):
        alias, separator, rest = node.partition("::")
        return mapping[alias] + separator + rest
    collection["modules"] = {mapping[alias]: value for alias, value in collection["modules"].items()}
    collection["module_paths"] = {mapping[alias]: value for alias, value in collection["module_paths"].items()}
    for key in ("items", "selected"):
        collection[key] = [rename(node) for node in collection[key]]
    collection["outcomes"] = {rename(node): value for node, value in collection["outcomes"].items()}
    tree = ET.parse(xml)
    for prop in tree.getroot().findall(".//property"):
        prop.set("value", rename(prop.get("value")))
    tree.write(xml)
    assert PROOF.proof_from_results(*checked_results)[0]["installed_wheel_tests"] == "passed"


@pytest.mark.parametrize("change", ["unknown", "duplicate_path", "missing_path", "selector", "config_hash", "config_path", "extra_flag"])
def test_actual_collection_or_explicit_config_mismatch_fails(checked_results, tmp_path, change):
    collection = checked_results[2]["collection"]
    aliases = list(collection["module_paths"])
    if change == "unknown":
        collection["module_paths"][aliases[0]] = str(tmp_path / "test_unknown.py")
    elif change == "duplicate_path":
        collection["module_paths"][aliases[1]] = collection["module_paths"][aliases[0]]
    elif change == "missing_path":
        collection.pop("module_paths")
    elif change == "selector":
        collection["selectors"]["python_functions"] = ["test_*"]
    elif change == "config_hash":
        collection["config"]["config_sha256"] = "0" * 64
    elif change == "config_path":
        copied = tmp_path / "different-pyproject.toml"
        copied.write_bytes(Path(collection["config"]["config_path"]).read_bytes())
        collection["config"]["config_path"] = str(copied)
    else:
        collection["invocation_args"].append("--ignore=some_case.py")
    with pytest.raises(ValueError):
        PROOF.proof_from_results(*checked_results)
