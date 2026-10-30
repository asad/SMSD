#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Collect a CUDA wheel only after retained XML results and installed identity agree."""

import argparse
import ast
import ctypes
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET
import zipfile

NATIVE_TESTS = {
    "smsd_batch_gpu_tests", "smsd_cuda_host_link_tests", "smsd_native_gpu_domain_tests",
    "smsd_gpu_domain_regression_tests", "smsd_screening_regression_tests", "smsd_screen_accelerator_tests",
}
_COLLECTION = {"items": set(), "modules": {}, "selected": [], "deselected": [], "outcomes": {}, "module_paths": {}}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Duplicate identity key: " + key)
            result[key] = value
        return result

    result = json.loads(Path(path).read_text(), object_pairs_hook=unique)
    if not isinstance(result, dict):
        raise ValueError("Expected a JSON object")
    return result


def junit_results(path):
    """Read actual cases, reject failures and retain every skipped identity."""
    root = ET.parse(path).getroot()
    suites = list(root.iter("testsuite"))
    cases = list(root.iter("testcase"))
    if not suites or not cases:
        raise ValueError("JUnit results contain no test suites or cases")
    for suite in ([root] + suites if root.tag == "testsuites" else suites):
        if suite.tag == "testsuite" and not suite.get("name"):
            raise ValueError("JUnit suite has no identity")
        for key in ("tests", "failures", "errors", "skipped"):
            value = suite.get(key, "0")
            if not re.fullmatch(r"[0-9]+", value):
                raise ValueError("Invalid JUnit count: " + key)
            if key in ("failures", "errors") and int(value):
                raise ValueError("JUnit reports " + key)
        nested_cases = list(suite.iter("testcase"))
        if (suite.tag == "testsuite" or "tests" in suite.attrib) and int(suite.get("tests", "0")) != len(nested_cases):
            raise ValueError("JUnit test count differs from its cases")
        skipped_count = sum(case.find("skipped") is not None or case.get("status") in ("notrun", "disabled")
                            for case in nested_cases)
        if "skipped" in suite.attrib and int(suite.get("skipped", "0")) != skipped_count:
            raise ValueError("JUnit skip count differs from recorded skipped cases")
    passed, skipped, identities, joined = [], [], set(), []
    for case in cases:
        name, classname = case.get("name", ""), case.get("classname", "")
        identity = classname + "::" + name
        if not name or identity in identities:
            raise ValueError("Missing or duplicate JUnit case identity")
        identities.add(identity)
        if (case.find("failure") is not None or case.find("error") is not None
                or case.get("status") in ("failed", "error")):
            raise ValueError("JUnit case failed: " + identity)
        skip = case.find("skipped")
        is_skipped = skip is not None or case.get("status") in ("notrun", "disabled")
        if is_skipped:
            skipped.append({"classname": classname, "name": name,
                            "reason": skip.get("message", "") if skip is not None else case.get("status")})
        elif not classname:
            raise ValueError("Passed JUnit case has no class identity")
        else:
            passed.append((classname, name))
        properties = [entry.get("value") for entry in case.findall("./properties/property")
                      if entry.get("name") == "smsd_nodeid"]
        if len(properties) > 1:
            raise ValueError("Duplicate XML node identity")
        joined.append({"classname": classname, "name": name, "nodeid": properties[0] if properties else None,
                       "outcome": "skipped" if is_skipped else "passed"})
    if not passed:
        raise ValueError("JUnit results contain no passed cases")
    return {"passed": passed, "skipped": skipped, "cases": joined}


def installed_files(wheel, package, prefix, module, native):
    """Join all installed package Python/native files to their exact wheel members."""
    prefix, module, native = Path(prefix).resolve(), Path(module).resolve(), Path(native).resolve()
    if not module.is_relative_to(prefix) or not native.is_relative_to(prefix):
        raise ValueError("Package was not loaded from the installed environment")
    base = module.parent.parent
    files = {}
    with zipfile.ZipFile(wheel) as archive:
        for member in archive.namelist():
            if (member.startswith(package + "/") or member.startswith(package + ".libs/")) and member.endswith((".py", ".so")):
                installed = (base / member).resolve()
                expected = hashlib.sha256(archive.read(member)).hexdigest()
                if not installed.is_relative_to(prefix) or sha256(installed) != expected:
                    raise ValueError("Installed package file differs from its wheel: " + member)
                files[member] = expected
    if (module != base / package / "__init__.py" or str(native.relative_to(base)) not in files
            or not files or {path.relative_to(base).as_posix() for path in module.parent.rglob("*.py")}
            != {name for name in files if name.startswith(package + "/") and name.endswith(".py")}):
        raise ValueError("Installed package origin or Python inventory differs from its wheel")
    return files


def pytest_collectstart(collector):
    if os.environ.get("SMSD_CUDA_TEST_IDENTITY") and collector.nodeid.endswith(".py"):
        path = str(Path(collector.path).resolve())
        previous = _COLLECTION["module_paths"].setdefault(collector.nodeid, path)
        if previous != path:
            raise ValueError("One collected module alias refers to different files")


def pytest_collectreport(report):
    if not os.environ.get("SMSD_CUDA_TEST_IDENTITY"):
        return
    if report.nodeid.endswith(".py"):
        _COLLECTION["modules"][report.nodeid] = {"outcome": report.outcome, "reason": str(report.longrepr)}
    for item in report.result or []:
        if callable(getattr(item, "runtest", None)):
            _COLLECTION["items"].add(item.nodeid)


def pytest_deselected(items):
    if os.environ.get("SMSD_CUDA_TEST_IDENTITY"):
        _COLLECTION["deselected"].extend(item.nodeid for item in items)


def pytest_collection_finish(session):
    if os.environ.get("SMSD_CUDA_TEST_IDENTITY"):
        _COLLECTION["selected"] = [item.nodeid for item in session.items]
        for item in session.items:
            item.user_properties.append(("smsd_nodeid", item.nodeid))


def pytest_runtest_logreport(report):
    if not os.environ.get("SMSD_CUDA_TEST_IDENTITY"):
        return
    if report.failed or report.skipped or report.when == "call":
        previous = _COLLECTION["outcomes"].get(report.nodeid, {})
        if previous.get("outcome") != "failed":
            _COLLECTION["outcomes"][report.nodeid] = {"outcome": report.outcome,
                                                       "reason": str(report.longrepr) if report.skipped else ""}


def expected_pytest_settings(config):
    """Use source-defined selectors and ordinary pytest defaults, not narrower globs."""
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib
    options = tomllib.loads(Path(config).read_text())["tool"]["pytest"]["ini_options"]
    def words(value):
        return shlex.split(value) if isinstance(value, str) else list(value)
    selectors = {name: words(options.get(name, default)) for name, default in (
        ("python_files", ["test_*.py", "*_test.py"]),
        ("python_classes", ["Test"]), ("python_functions", ["test"]))}
    return selectors, words(options.get("addopts", "")), str(options.get("timeout", ""))


def module_paths(collection, test_roots):
    """Join actual collector paths bijectively to the complete declared suite files."""
    expected = {path.resolve() for directory in test_roots for pattern in ("test_*.py", "*_test.py")
                for path in directory.rglob(pattern)}
    recorded = collection.get("module_paths")
    if recorded is None or set(recorded) != set(collection["modules"]):
        raise ValueError("Collected module paths are missing or have unknown aliases")
    joined = {}
    for alias, filename in recorded.items():
        path = Path(filename).resolve()
        if (not alias.endswith(".py") or "::" in alias or path not in expected
                or Path(alias).name != path.name):
            raise ValueError("Collected module path or alias is unknown")
        joined[alias] = path
    if not expected or len(set(joined.values())) != len(joined) or set(joined.values()) != expected:
        raise ValueError("Python module collection differs from the complete source suites")
    return joined


def loaded_capabilities(cuda):
    """Read capability metadata from the same installed library, without a screen."""
    available = cuda.is_available()
    function = cuda._library.smsd_cuda_binary_mma_available
    function.argtypes, function.restype = [ctypes.c_void_p, ctypes.c_size_t], ctypes.c_int
    error = ctypes.create_string_buffer(256)
    mma = function(error, len(error))
    if mma < 0:
        raise ValueError("Could not query binary MMA capability: " + error.value.decode())
    return {"available": available, "device_info": cuda.device_info(),
            "binary_mma_available": mma, "binary_mma_error": error.value.decode(),
            "testing_failure_hook_present": hasattr(cuda._library, "smsd_cuda_test_fail_next")}


def checked_capabilities(identity):
    capabilities = identity.get("capabilities")
    if (not capabilities or capabilities.get("available") is not True
            or not re.search(r"compute [0-9]+\.[0-9]+", capabilities.get("device_info", ""))
            or type(capabilities.get("binary_mma_available")) is not int
            or capabilities["binary_mma_available"] not in (0, 1)
            or capabilities.get("binary_mma_error") != ""
            or type(capabilities.get("testing_failure_hook_present")) is not bool):
        raise ValueError("Positive installed-device capability evidence is missing")
    return capabilities


def capability_skip(name, reason, actual_reason, capabilities):
    """Only named, parameter-exact unsupported-feature controls may be skipped."""
    mma = {"test_opt_in_binary_mma_integer_intersections[" + str(words) + "]"
           for words in (1, 3, 4, 5, 33)}
    injection = {function + "[" + str(failure) + "]" for failure in (1, 2) for function in (
        "test_resident_top_k_injected_failure_keeps_payload_and_handle",
        "test_resident_graph_failure_recovery_preserves_payload")}
    expected = None
    if name in mma and capabilities["binary_mma_available"] == 0:
        expected = "binary MMA prototype needs an enabled build and an SM80-SM90 device"
    elif (name == "test_binary_mma_request_keeps_unsupported_device_on_portable_path"
          and capabilities["binary_mma_available"] > 0):
        expected = "this control requires the portable fallback device or build"
    elif name in injection and not capabilities["testing_failure_hook_present"]:
        expected = "failure injection requires the separate SMSD_CUDA_TESTING library"
    try:
        parsed = ast.literal_eval(actual_reason)
    except (ValueError, SyntaxError):
        parsed = None
    if isinstance(parsed, tuple) and len(parsed) == 3 and isinstance(parsed[2], str):
        actual_reason = parsed[2].removeprefix("Skipped: ")
    if expected is None or reason != expected or actual_reason != expected:
        raise ValueError("Required CUDA case was skipped or its capability reason differs")


def check_collection(identity, python):
    collection = identity["collection"]
    if (collection["deselected"] or collection["filters"] or collection["env_addopts"]
            or len(collection["selected"]) != len(set(collection["selected"]))
            or set(collection["items"]) != set(collection["selected"])
            or set(collection["outcomes"]) != set(collection["selected"])):
        raise ValueError("Python collection was filtered, deselected or truncated")
    root = Path(collection["rootdir"])
    test_roots = [Path(path).resolve() for path in identity["test_roots"]]
    config_path = test_roots[0].parents[1] / "pyproject.toml"
    selectors, addopts, timeout = expected_pytest_settings(config_path)
    expected_args = [str(path) for path in test_roots]
    expected_invocation = [*expected_args, "-q", "-rs", "-p", "no:cacheprovider", "-p", "cuda_wheel_provenance",
                           "--junitxml=" + identity["junit_xml"], "-c", str(config_path), "--rootdir", str(root)]
    config = collection.get("config")
    if collection["invocation_args"] != expected_invocation or config is None:
        raise ValueError("Python invocation does not select the canonical pytest configuration")
    if (collection["args"] != expected_args or collection["selectors"] != selectors
            or collection["ini_addopts"] != addopts or addopts != ["--import-mode=importlib"]
            or timeout != "60" or config["config_sha256"] != sha256(config_path)
            or sha256(config["config_path"]) != sha256(config_path)
            or config["importmode"] != "importlib" or str(config["timeout"]) != timeout
            or config["ini_addopts"] != addopts or config["invocation_args"] != collection["invocation_args"]
            or Path(config["rootdir"]).resolve() != root.resolve()):
        raise ValueError("Python configuration differs from the complete canonical suites")
    if Path(config["config_path"]).resolve() != config_path.resolve():
        raise ValueError("Explicit pytest configuration differs from the SDK")
    paths = module_paths(collection, test_roots)
    capabilities = checked_capabilities(identity)
    collection_skips = {(alias.replace("/", ".").removesuffix(".py"),)
                        for alias, report in collection["modules"].items() if report["outcome"] == "skipped"}
    for alias, report in collection["modules"].items():
        if report["outcome"] not in ("passed", "skipped") or (report["outcome"] == "skipped" and not report["reason"]):
            raise ValueError("Python module did not collect successfully")
        if paths[alias].is_relative_to(test_roots[1]) and report["outcome"] != "passed":
            raise ValueError("Required CUDA module was skipped")
    skips = {(entry["classname"], entry["name"]): entry["reason"] for entry in python["skipped"]}
    joined, seen_collection_skips = {}, set()
    for case in python["cases"]:
        nodeid = case["nodeid"]
        if nodeid is None:
            key = (case["name"],)
            if case["outcome"] != "skipped" or case["classname"] or key not in collection_skips:
                raise ValueError("XML case has no matching collected node")
            seen_collection_skips.add(key)
            continue
        alias = nodeid.split("::", 1)[0]
        if alias not in paths or nodeid in joined or nodeid not in collection["outcomes"]:
            raise ValueError("XML node is duplicate or absent from collection")
        actual = collection["outcomes"][nodeid]
        if actual["outcome"] != case["outcome"] or actual["outcome"] not in ("passed", "skipped"):
            raise ValueError("XML outcome differs from the actual test report")
        if actual["outcome"] == "skipped" and not actual["reason"]:
            raise ValueError("A default skip has no retained context")
        if paths[alias].is_relative_to(test_roots[1]) and case["outcome"] == "skipped":
            if paths[alias].name != "test_resident_fingerprints.py":
                raise ValueError("Required CUDA case was skipped")
            capability_skip(case["name"], skips[(case["classname"], case["name"])], actual["reason"], capabilities)
        joined[nodeid] = case["outcome"]
    if set(joined) != set(collection["selected"]) or seen_collection_skips != collection_skips:
        raise ValueError("XML results are missing collected cases or module skips")


def loaded_runtime(cuda):
    """Retain the runtime actually loaded by the CUDA tests, including its library bytes."""
    import importlib.metadata

    runtime_version = ctypes.c_int()
    get_version = cuda._library.cudaRuntimeGetVersion
    get_version.argtypes = [ctypes.POINTER(ctypes.c_int)]
    get_version.restype = ctypes.c_int
    if get_version(ctypes.byref(runtime_version)) != 0:
        raise ValueError("Could not read the loaded CUDA runtime version")
    paths = {Path(fields[5]).resolve() for line in Path("/proc/self/maps").read_text().splitlines()
             if len(fields := line.split(maxsplit=5)) == 6 and "libcudart.so" in fields[5]}
    if len(paths) != 1:
        raise ValueError("The loaded CUDA runtime library is missing or ambiguous")
    return {"version": runtime_version.value,
            "package_version": importlib.metadata.version("nvidia-cuda-runtime-cu12"),
            "libraries": {str(path): sha256(path) for path in paths}}


def check_installed_identity(identity, source, wheel):
    """Check the session-end files against the exact repaired wheel bytes."""
    if (type(identity.get("pytest_exit_status")) is not int or identity["pytest_exit_status"] != 0
            or identity.get("source_commit") != source["source_commit"]
            or identity.get("sdist_sha256") != source["sdist_sha256"]
            or identity.get("smsd_version") != source["version"]
            or identity.get("wheel") != wheel.name or identity.get("wheel_sha256") != sha256(wheel)):
        raise ValueError("Test identity differs from the successful source/wheel invocation")
    prefix = Path(identity["prefix"]).resolve()
    with zipfile.ZipFile(wheel) as archive:
        for field, member in (("module", "smsd_cuda/__init__.py"),
                              ("library", "smsd_cuda/libsmsd_cuda_screen.so")):
            installed = Path(identity[field]).resolve()
            if not installed.is_relative_to(prefix):
                raise ValueError("CUDA module was not loaded from the installed environment")
            expected = hashlib.sha256(archive.read(member)).hexdigest()
            if identity.get(field + "_sha256") != expected or sha256(installed) != expected:
                raise ValueError("Installed CUDA file differs from the tested repaired wheel")
    native = Path(identity["native_module"]).resolve()
    if not native.is_relative_to(prefix) or sha256(native) != identity["native_module_sha256"]:
        raise ValueError("Installed core native module identity differs")
    cpu_wheel = Path(identity["cpu_wheel"])
    if sha256(cpu_wheel) != identity["cpu_wheel_sha256"]:
        raise ValueError("CPU wheel changed after its test invocation")
    if installed_files(cpu_wheel, "smsd", prefix, identity["core_module"], native) != identity["core_files"]:
        raise ValueError("Installed core files differ from the exact CPU wheel")
    if installed_files(wheel, "smsd_cuda", prefix, identity["module"], identity["library"]) != identity["cuda_files"]:
        raise ValueError("Installed CUDA files differ from the exact CUDA wheel")


def pytest_sessionfinish(session, exitstatus):
    """Retain identity from the same ordinary pytest process that loaded the packages."""
    output = os.environ.get("SMSD_CUDA_TEST_IDENTITY")
    if not output:
        return
    wheel = Path(os.environ["SMSD_CUDA_TEST_WHEEL"]).resolve()
    manifest = Path(os.environ["SMSD_CUDA_SOURCE_MANIFEST"]).resolve()
    source = read_json(manifest)
    build_note = read_json(Path(os.environ["SMSD_CUDA_BUILD_NOTE"]))
    core, cuda = sys.modules.get("smsd"), sys.modules.get("smsd_cuda")
    if core is None or cuda is None or getattr(cuda, "_library", None) is None:
        raise ValueError("Tests did not load the installed core and CUDA library")
    if (core.__version__ != source["version"] or core._smsd.__version__ != source["version"]
            or cuda.__version__ != source["version"] or cuda.__smsd_version__ != source["version"]):
        raise ValueError("Loaded package versions differ from the canonical source")
    module = Path(cuda.__file__).resolve()
    library = Path(cuda._library._name).resolve()
    native = Path(core._smsd.__file__).resolve()
    filters = {name: getattr(session.config.option, name, None)
               for name in ("keyword", "markexpr", "deselect", "ignore", "ignore_glob", "lf", "ff", "stepwise",
                            "collectonly", "maxfail", "pyargs", "continue_on_collection_errors")
               if getattr(session.config.option, name, None)}
    collection = {**_COLLECTION, "items": sorted(_COLLECTION["items"]),
                  "rootdir": str(session.config.rootpath), "filters": filters,
                  "env_addopts": os.environ.get("PYTEST_ADDOPTS", ""),
                  "ini_addopts": session.config.getini("addopts"),
                  "invocation_args": list(session.config.invocation_params.args),
                  "args": [str(argument) for argument in session.config.args],
                  "config": {"config_path": str(session.config.inipath.resolve()),
                             "config_sha256": sha256(session.config.inipath),
                             "importmode": session.config.option.importmode,
                             "timeout": str(session.config.getini("timeout")),
                             "ini_addopts": session.config.getini("addopts"),
                             "invocation_args": list(session.config.invocation_params.args),
                             "rootdir": str(session.config.rootpath)},
                  "pytest_version": __import__("importlib.metadata", fromlist=["version"]).version("pytest"),
                  "pytest_config_source": {"path": str(Path(sys.modules["_pytest.config"].__file__).resolve()),
                                           "sha256": sha256(sys.modules["_pytest.config"].__file__)},
                  "selectors": {name: session.config.getini(name)
                                for name in ("python_files", "python_classes", "python_functions")}}
    identity = {
        "source_commit": source["source_commit"],
        "sdist_sha256": sha256(manifest.parent / source["sdist"]),
        "smsd_version": source["version"], "wheel": wheel.name, "wheel_sha256": sha256(wheel),
        "prefix": str(Path(sys.prefix).resolve()), "module": str(module), "module_sha256": sha256(module),
        "library": str(library), "library_sha256": sha256(library),
        "native_module": str(native), "native_module_sha256": sha256(native),
        "pytest_exit_status": int(exitstatus), "junit_xml": str(Path(session.config.option.xmlpath).resolve()),
        "core_module": str(Path(core.__file__).resolve()), "cpu_wheel": build_note["cpu_wheel"],
        "cpu_wheel_sha256": sha256(build_note["cpu_wheel"]),
        "core_files": installed_files(Path(build_note["cpu_wheel"]), "smsd", sys.prefix, core.__file__, native),
        "cuda_files": installed_files(wheel, "smsd_cuda", sys.prefix, module, library),
        "build_note": str(Path(os.environ["SMSD_CUDA_BUILD_NOTE"]).resolve()),
        "build_note_sha256": sha256(os.environ["SMSD_CUDA_BUILD_NOTE"]),
        "native_xml": str(Path(os.environ["SMSD_CUDA_NATIVE_XML"]).resolve()),
        "native_xml_sha256": sha256(os.environ["SMSD_CUDA_NATIVE_XML"]),
        "ctest_exit_status": int(os.environ["SMSD_CUDA_CTEST_STATUS"]),
        "ctest_started_ns": int(os.environ["SMSD_CUDA_CTEST_STARTED_NS"]),
        "ctest_finished_ns": int(os.environ["SMSD_CUDA_CTEST_FINISHED_NS"]),
        "runtime": loaded_runtime(cuda), "capabilities": loaded_capabilities(cuda),
        "collection": collection, "test_roots": [str(Path(build_note["sdk_source"]) / "python/tests"),
                                                   str(Path(build_note["checkout"]) / "packaging/smsd-cuda/tests")],
    }
    with Path(output).open("x") as stream:
        json.dump(identity, stream, indent=2)
        stream.write("\n")


def proof_from_results(source, wheel, identity, native_xml, python_xml):
    check_installed_identity(identity, source, wheel)
    if Path(identity["junit_xml"]).resolve() != python_xml.resolve():
        raise ValueError("Python XML differs from the test invocation")
    note_path = Path(identity["build_note"])
    note = read_json(note_path)
    if (sha256(note_path) != identity["build_note_sha256"] or note["source_commit"] != source["source_commit"]
            or note["sdist_sha256"] != source["sdist_sha256"]
            or note["cpu_wheel_sha256"] != identity["cpu_wheel_sha256"]
            or note["cuda_wheel_sha256"] != identity["wheel_sha256"]
            or type(identity["ctest_exit_status"]) is not int or identity["ctest_exit_status"] != 0
            or Path(identity["native_xml"]).resolve() != native_xml.resolve()
            or sha256(native_xml) != identity["native_xml_sha256"]
            or not note["finished_ns"] <= identity["ctest_started_ns"] <= native_xml.stat().st_mtime_ns
                   <= identity["ctest_finished_ns"]):
        raise ValueError("Native XML or build identity is stale or mismatched")
    if (identity["test_roots"] != [str(Path(note["sdk_source"]) / "python/tests"),
                                   str(Path(note["checkout"]) / "packaging/smsd-cuda/tests")]
            or Path(identity["collection"]["rootdir"]).resolve() != Path(note["sdk_source"]).resolve()):
        raise ValueError("Python suite roots differ from the checked build inputs")
    if set(note["native_binaries"]) != NATIVE_TESTS or any(
            sha256(Path(note["native_dir"]) / name) != digest for name, digest in note["native_binaries"].items()):
        raise ValueError("Native binaries differ from the checked build")
    runtime = identity["runtime"]
    if (not re.search(r"release 12\.8,", note["cuda_toolkit"]) or runtime["version"] != 12080
            or not re.match(r"12\.8(?:\.|$)", runtime["package_version"])
            or len(runtime["libraries"]) != 1
            or any(sha256(path) != digest for path, digest in runtime["libraries"].items())):
        raise ValueError("Actual toolkit/runtime evidence does not cover the declared CUDA 12.8 floor")
    native = junit_results(native_xml)
    if native["skipped"] or len(native["passed"]) != len(NATIVE_TESTS) or {name for _, name in native["passed"]} != NATIVE_TESTS:
        raise ValueError("Native XML does not contain the complete passed CUDA suite selection")
    python = junit_results(python_xml)
    check_collection(identity, python)
    if not any("test_smsd" in classname.split(".") for classname, _ in python["passed"]):
        raise ValueError("Python XML has no passed core suite")
    if not any("test_smsd_cuda" in classname.split(".") and name == "test_smsd_attaches_the_device"
               for classname, name in python["passed"]):
        raise ValueError("Python XML has no passed CUDA device test")
    if not any("test_resident_fingerprints" in classname.split(".")
               and name.startswith("test_resident_integer_thresholds_and_order[")
               for classname, name in python["passed"]):
        raise ValueError("Python XML has no passed resident CUDA callback test")
    proof = {"schema_version": 1, "source_commit": source["source_commit"],
             "sdist_sha256": source["sdist_sha256"], "smsd_version": source["version"],
             "wheel": wheel.name, "wheel_sha256": sha256(wheel),
             "installed_wheel_tests": "passed", "native_gpu_tests": "passed"}
    summary = {"native_passed": len(native["passed"]), "python_passed": len(python["passed"]),
               "python_skipped": python["skipped"], "native_xml_sha256": sha256(native_xml),
               "python_xml_sha256": sha256(python_xml), "collection": identity["collection"],
               "runtime": runtime, "cuda_toolkit": note["cuda_toolkit"],
               "capabilities": checked_capabilities(identity)}
    return proof, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--identity", type=Path, required=True)
    parser.add_argument("--native-xml", type=Path, required=True)
    parser.add_argument("--python-xml", type=Path, required=True)
    parser.add_argument("--evidence-summary", type=Path, required=True)
    parser.add_argument("--checkout", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    release_root = args.release_dir.resolve()
    proof_path = args.identity.with_name("cuda-wheel-build-proof.json")
    identity = read_json(args.identity)
    for private in (args.identity, args.native_xml, args.python_xml, args.evidence_summary,
                    proof_path, Path(identity["build_note"])):
        if private.resolve().is_relative_to(release_root):
            raise ValueError("Private validation evidence must be outside the release directory")
    spec = importlib.util.spec_from_file_location("cuda_collector", Path(__file__).with_name("collect_cuda_wheels.py"))
    collector = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(collector)
    manifest = args.release_dir / "SOURCE_MANIFEST.json"
    source = collector.BUILD.read_manifest(manifest)
    source = collector.BUILD.checked_source(args.release_dir / source["sdist"], manifest)
    note = read_json(identity["build_note"])
    if Path(note["checkout"]).resolve() != args.checkout.resolve() or Path(note["source_manifest"]).resolve() != manifest.resolve():
        raise ValueError("Build inputs differ from the source/checkout being collected")
    subprocess.run([sys.executable, str(args.checkout / "scripts/build_python_wheels.py"), "check-source",
                    "--sdist", str((args.release_dir / source["sdist"]).resolve()),
                    "--source-manifest", str(manifest.resolve()), "--release-commit", source["source_commit"]],
                   cwd=args.checkout, check=True, stdout=subprocess.DEVNULL)
    expected = collector.BUILD.read_source_archive(args.release_dir / source["sdist"])
    sdk = Path(note["sdk_source"])
    actual = {path.relative_to(sdk).as_posix() for path in sdk.rglob("*")
              if path.is_file() and "__pycache__" not in path.parts and path.suffix != ".pyc"}
    if actual != set(expected) or any((sdk / name).read_bytes() != data for name, data in expected.items()):
        raise ValueError("Test SDK differs from the canonical source archive")
    proof, summary = proof_from_results(source, args.wheel, identity, args.native_xml, args.python_xml)
    # Paths, case details and XML remain local; only the small checked proof enters the release.
    collector.BUILD.write_manifest(args.evidence_summary, summary)
    collector.BUILD.write_manifest(proof_path, proof)
    collector.collect(args.release_dir, args.wheel, proof_path, args.checkout)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, ET.ParseError, zipfile.BadZipFile, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error)) from error
