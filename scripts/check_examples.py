#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Run the standalone Python, Java and C++ examples in the current guides.

The Python checks also run every script and notebook in python/examples, and
drive each tool-calling loop in docs/AI_ASSISTANTS.md with fixed model replies,
so that no network, model service or API key is needed.
"""

import argparse
import ast
import contextlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import types


ROOT = Path(__file__).resolve().parents[1]
GUIDES = (
    "README.md",
    "python/README.md",
    "java/README.md",
    "cpp/README.md",
    "docs/PYTHON.md",
    "docs/JAVA.md",
    "docs/CPP.md",
    "docs/EXAMPLES.md",
    "docs/PACKED_CORPORA.md",
    "docs/SMARTS_EXPORT.md",
    "docs/PYTHON_ADAPTERS.md",
    "docs/guide/README.md",
    "docs/guide/GETTING_STARTED.md",
    "docs/guide/SUBSTRUCTURE.md",
    "docs/guide/SMARTS.md",
    "docs/guide/MCS.md",
    "docs/guide/FINGERPRINTS.md",
    "docs/guide/SCREENING.md",
    "docs/guide/PACKED_INDEX.md",
    "docs/guide/CPU_GPU.md",
    "docs/guide/COMMAND_LINE.md",
    "docs/guide/OTHER_LANGUAGES.md",
)
ASSISTANT_GUIDE = "docs/AI_ASSISTANTS.md"
EXAMPLE_DIRECTORY = "python/examples"
# Optional packages named by the smsd extras and the example scripts; a script whose
# run stops on one of these, when absent, is skipped rather than failed.
OPTIONAL_MODULES = ("IPython", "mcp", "numpy", "pandas", "polars", "pyarrow", "rdkit", "scipy",
                    "sklearn", "smsd_cuda", "torch", "torch_geometric")
FENCE = re.compile(r"^```(python|java|cpp)[ \t]*\n(.*?)^```[ \t]*$", re.M | re.S)
IMPORT_FAILURE = re.compile(r"^(?:\w+\.)*(?:ModuleNotFoundError|ImportError): (.*)$")

# The fixed model turn asks for one valid and one invalid call; the next turn ends the loop.
TOOL_CALLS = (("call_1", "find_mcs", {"query_smiles": "c1ccccc1O", "target_smiles": "c1ccccc1N"}),
              ("call_2", "find_mcs", {"query_smiles": "c1ccccc1O"}))
FINAL_TEXT = "Offline tool loop finished."
# Model endpoints by the last two attribute names of the call, and the reply style each uses.
ENDPOINTS = {("messages", "create"): "messages", ("models", "generate_content"): "contents",
             ("responses", "create"): "responses", ("completions", "create"): "chat"}
HARNESS = ("import importlib.util, sys; "
           "spec = importlib.util.spec_from_file_location('check_examples', sys.argv[1]); "
           "module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); "
           "module.offline_tool_loop(sys.argv[2])")


def run(command, directory, timeout, extra_environment=None):
    environment = os.environ.copy()
    environment.pop("PYTHONOPTIMIZE", None)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment.update(extra_environment or {})
    result = subprocess.run(
        [str(part) for part in command],
        cwd=directory,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        env=environment,
    )
    if result.returncode:
        output = (result.stdout + result.stderr).strip()
        raise RuntimeError(f"Command exited {result.returncode}: {output}")
    return result.stdout.strip()


def examples(languages):
    for guide in GUIDES:
        text = (ROOT / guide).read_text(encoding="utf-8")
        for match in FENCE.finditer(text):
            language, code = match.groups()
            if language in languages:
                yield guide, text.count("\n", 0, match.start()) + 1, language, code


def third_party(name):
    """The top-level module of an import outside smsd and the standard library, or None."""
    top = name.partition(".")[0]
    return None if top == "smsd" or top in sys.stdlib_module_names else top


def imported(node):
    if isinstance(node, ast.Import):
        return [alias.name for alias in node.names]
    if isinstance(node, ast.ImportFrom) and not node.level and node.module:
        return [node.module]
    return []


def falls_back(statement):
    """Whether a try statement carries on without the module when its import fails."""
    for handler in statement.handlers:
        caught = [] if handler.type is None else getattr(handler.type, "elts", [handler.type])
        if handler.type is not None and not any(
                isinstance(name, ast.Name) and name.id in ("ImportError", "ModuleNotFoundError", "Exception")
                for name in caught):
            continue
        return not any(isinstance(node, ast.Raise) or isinstance(node, ast.Call) and (
            isinstance(node.func, ast.Attribute) and node.func.attr == "exit"
            or isinstance(node.func, ast.Name) and node.func.id in ("exit", "quit", "SystemExit"))
            for node in ast.walk(handler))
    return False


def required_modules(source):
    """Third-party modules a script needs before it can start.

    These are the module-level imports, except those in a try statement that
    carries on without them. Imports under if, def or class are judged when run.
    """
    required = set()

    def visit(statements):
        for statement in statements:
            if isinstance(statement, ast.Try):
                if not falls_back(statement):
                    visit(statement.body)
            else:
                required.update(filter(None, map(third_party, imported(statement))))

    visit(ast.parse(source).body)
    return required


def missing_optional(error, absent):
    """The absent optional module named by the final import error of a failed run."""
    lines = [line.strip() for line in str(error).splitlines() if line.strip()]
    match = IMPORT_FAILURE.match(lines[-1]) if lines else None
    if match:
        for name in sorted(absent):
            if re.search(r"(?<![\w.])" + re.escape(name) + r"(?!\w)", match.group(1), re.I):
                return name
    return None


def notebook_source(path):
    """The code cells of a notebook as one script, or None if a cell uses IPython syntax."""
    parts = []
    for number, cell in enumerate(json.loads(path.read_text(encoding="utf-8"))["cells"], 1):
        if cell.get("cell_type") != "code":
            continue
        source = cell["source"] if isinstance(cell["source"], str) else "".join(cell["source"])
        if re.search(r"^[ \t]*[%!]", source, re.M):
            return None
        parts.append(f"# Cell {number}\n{source}\n")
    return "\n".join(parts)


def example_scripts():
    """Each example script or notebook with its Python source."""
    directory = ROOT / EXAMPLE_DIRECTORY
    for path in sorted(directory.glob("*.py")):
        yield path, path.read_text(encoding="utf-8")
    for path in sorted(directory.glob("*.ipynb")):
        yield path, notebook_source(path)


def assistant_loops():
    path = ROOT / ASSISTANT_GUIDE
    text = path.read_text(encoding="utf-8") if path.is_file() else ""
    for match in FENCE.finditer(text):
        language, code = match.groups()
        if language == "python":
            yield ASSISTANT_GUIDE, text.count("\n", 0, match.start()) + 1, code


class FakeName:
    """Any name from a model SDK: attributes are further names, and calling one returns
    an object holding the keyword arguments, except at a model endpoint, which answers
    from the fixed replies."""

    def __init__(self, path, model, values=None):
        self._path, self._model = path, model
        self.__dict__.update(values or {})

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        return FakeName(self._path + (name,), self._model)

    def __call__(self, *args, **kwargs):
        style = ENDPOINTS.get(self._path[-2:])
        if style:
            return self._model.reply(style, kwargs)
        if self._path[-1] in ("create", "generate_content", "stream"):
            raise RuntimeError("No offline reply for " + ".".join(self._path))
        return FakeName(self._path, self._model, kwargs)


class FakeModel:
    """Asks for the fixed tool calls, then ends the conversation."""

    def __init__(self):
        self.requests = []

    def reply(self, style, request):
        self.requests.append((style, {key: list(value) if isinstance(value, list) else value
                                      for key, value in request.items()}))
        if len(self.requests) > 2:
            raise RuntimeError("The loop asked the model again after the final reply")
        Reply = types.SimpleNamespace
        if len(self.requests) == 2:
            return {"messages": Reply(stop_reason="end_turn", content=[Reply(type="text", text=FINAL_TEXT)]),
                    "contents": Reply(function_calls=None, text=FINAL_TEXT,
                                      candidates=[Reply(content=Reply(role="model", parts=[]))]),
                    "responses": Reply(output=[Reply(type="message", content=[])], output_text=FINAL_TEXT),
                    "chat": Reply(choices=[Reply(message=Reply(role="assistant", content=FINAL_TEXT,
                                                               tool_calls=None))])}[style]
        if style == "messages":
            return Reply(stop_reason="tool_use", content=[Reply(type="text", text="Searching.")] + [
                Reply(type="tool_use", id=key, name=name, input=dict(arguments))
                for key, name, arguments in TOOL_CALLS])
        if style == "contents":
            return Reply(text=None, candidates=[Reply(content=Reply(role="model", parts=[]))],
                         function_calls=[Reply(id=key, name=name, args=dict(arguments))
                                         for key, name, arguments in TOOL_CALLS])
        if style == "responses":
            return Reply(output_text="", output=[
                Reply(type="function_call", call_id=key, name=name, arguments=json.dumps(arguments))
                for key, name, arguments in TOOL_CALLS])
        return Reply(choices=[Reply(message=Reply(role="assistant", content=None, tool_calls=[
            Reply(id=key, type="function", function=Reply(name=name, arguments=json.dumps(arguments)))
            for key, name, arguments in TOOL_CALLS]))])


def returned_results(style, request):
    """(call id, result, error) for each tool result the loop sent back to the model."""
    def decoded(payload, flagged=None):
        if isinstance(payload, dict):
            return payload.get("result"), payload.get("error")
        try:
            value = json.loads(payload)
        except (TypeError, ValueError):
            value = None
        if isinstance(value, dict) and flagged is not True:
            return value, None
        return None, str(payload)

    if style == "messages":
        return [(item["tool_use_id"], *decoded(item["content"], item.get("is_error")))
                for item in request["messages"][-1]["content"]]
    if style == "contents":
        return [(part.function_response.id, *decoded(part.function_response.response))
                for part in request["contents"][-1].parts]
    if style == "responses":
        return [(item["call_id"], *decoded(item["output"])) for item in request["input"]
                if isinstance(item, dict) and item.get("type") == "function_call_output"]
    return [(item["tool_call_id"], *decoded(item["content"])) for item in request["messages"]
            if isinstance(item, dict) and item.get("role") == "tool"]


def mentions(value, text, depth=0):
    if depth > 12:
        return False
    if isinstance(value, str):
        return value == text
    if isinstance(value, FakeName):
        value = {key: item for key, item in vars(value).items() if not key.startswith("_")}
    if isinstance(value, dict):
        return any(mentions(item, text, depth + 1) for item in value.values())
    if isinstance(value, (list, tuple)):
        return any(mentions(item, text, depth + 1) for item in value)
    return False


def fake_module(dotted, model):
    """A stand-in package whose every public name is a FakeName."""
    module = types.ModuleType(dotted)

    def attribute(name):
        if name.startswith("__"):
            raise AttributeError(name)
        return FakeName((dotted, name), model)

    module.__getattr__ = attribute
    module.__path__, module.offline_fake = [], True
    return module


def offline_tool_loop(path):
    """Run one published tool loop with fixed model replies and check what it sent back."""
    import socket

    def refuse(*args, **kwargs):
        raise OSError("Network access is disabled in the offline tool-loop check")

    socket.socket.connect = socket.socket.connect_ex = refuse
    socket.create_connection = socket.getaddrinfo = refuse
    from smsd.tool_definitions import call_tool

    code = Path(path).read_text(encoding="utf-8")
    model = FakeModel()
    names = {name for node in ast.walk(ast.parse(code)) for name in imported(node) if third_party(name)}
    if not names:
        raise SystemExit("The loop imports no model client to replace")
    for name in sorted(names):
        parts = name.split(".")
        for end in range(1, len(parts) + 1):
            dotted = ".".join(parts[:end])
            if not getattr(sys.modules.get(dotted), "offline_fake", False):
                sys.modules[dotted] = fake_module(dotted, model)
            if end > 1:
                setattr(sys.modules[".".join(parts[:end - 1])], parts[end - 1], sys.modules[dotted])

    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        exec(compile(code, path, "exec"), {"__name__": "__main__"})
    if len(model.requests) != 2:
        raise SystemExit(f"Expected two model requests, saw {len(model.requests)}")
    style, first = model.requests[0]
    if not mentions(first, TOOL_CALLS[0][1]):
        raise SystemExit("The first request did not offer the SMSD tools")
    results = returned_results(style, model.requests[1][1])
    if [result[0] for result in results] != [key for key, _, _ in TOOL_CALLS]:
        raise SystemExit(f"Tool results do not answer the calls in order: {results!r}")
    (_, value, error), (_, invalid, message) = results
    expected = call_tool(TOOL_CALLS[0][1], dict(TOOL_CALLS[0][2]))
    fields = ("found", "status", "optimal", "complete", "score")
    if error or not isinstance(value, dict) or any(value.get(field) != expected.get(field) for field in fields):
        raise SystemExit(f"The valid call did not return the local tool result: {value!r} {error!r}")
    try:
        call_tool(TOOL_CALLS[1][1], dict(TOOL_CALLS[1][2]))
    except ValueError as reason:
        expected_error = str(reason)
    else:
        raise SystemExit("The invalid fixed call was accepted by call_tool")
    if invalid is not None or expected_error not in (message or ""):
        raise SystemExit(f"The invalid call did not return its error to the model: {message!r}")
    if FINAL_TEXT not in output.getvalue():
        raise SystemExit("The loop did not print the final model text")
    print(json.dumps({"reply_style": style, "tool_results": len(results), "fields_checked": fields}))


def default_jar():
    metadata = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    version = re.search(r'^version = "([^"]+)"', metadata, re.M).group(1)
    filename = f"smsd-cli-{version}-jar-with-dependencies.jar"
    for directory in (ROOT / "java/cli/target", ROOT / "dist" / f"release-{version}"):
        candidate = directory / filename
        if candidate.is_file():
            return candidate
    return None


def report(results, record):
    results.append(record)
    kind = "" if record["kind"] == "code block" else ", " + record["kind"]
    print(f"{record['status'].upper()}: {record['guide']}:{record['line']} ({record['language']}{kind})",
          flush=True)
    if record.get("reason"):
        print(f"  {record['reason']}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--languages", nargs="+", choices=("python", "java", "cpp"),
        default=["python", "java", "cpp"],
    )
    parser.add_argument("--jar", type=Path, help="SMSD JAR including its dependencies")
    parser.add_argument("--java-home", type=Path, help="JDK used to compile and run Java examples")
    parser.add_argument("--include-dir", type=Path, default=ROOT / "cpp/include")
    parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"), help="GCC or Clang C++ compiler command")
    parser.add_argument("--timeout", type=float, default=120, help="Seconds allowed per command")
    parser.add_argument("--report", type=Path, help="Write a JSON result report")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")

    jar = args.jar.resolve() if args.jar else default_jar()
    java_home = args.java_home or os.environ.get("JAVA_HOME")
    if java_home:
        java_home = Path(java_home).resolve()
    suffix = ".exe" if os.name == "nt" else ""
    java = Path(java_home) / "bin" / f"java{suffix}" if java_home else shutil.which("java")
    javac = Path(java_home) / "bin" / f"javac{suffix}" if java_home else shutil.which("javac")
    executable = args.cxx.strip('"')
    compiler = [executable] if Path(executable).is_file() else shlex.split(args.cxx)
    if compiler and shutil.which(compiler[0]):
        compiler[0] = os.path.abspath(shutil.which(compiler[0]))
    results = []
    environment = {}

    with tempfile.TemporaryDirectory(prefix="smsd-examples-") as temporary:
        base = Path(temporary)
        rdkit_available = False
        scripts, absent = [], set()
        if "python" in args.languages:
            scripts = list(example_scripts())
            candidates = set(OPTIONAL_MODULES)
            for _, source in scripts:
                try:
                    candidates.update(required_modules(source or ""))
                except SyntaxError:
                    pass
            try:
                probe = run(
                    [sys.executable, "-c", "import importlib.util, json, sys, smsd; "
                     "print(json.dumps({'python': sys.version.split()[0], "
                     "'smsd': smsd.__version__, 'smsd_path': smsd.__file__, "
                     "'rdkit_available': importlib.util.find_spec('rdkit') is not None, "
                     "'optional_modules': {name: importlib.util.find_spec(name) is not None "
                     "for name in sys.argv[1:]}}))", *sorted(candidates)],
                    base, args.timeout,
                )
                environment = json.loads(probe)
                rdkit_available = environment["rdkit_available"]
                absent = {name for name, found in environment["optional_modules"].items() if not found}
            except (RuntimeError, subprocess.TimeoutExpired) as error:
                parser.exit(1, f"Python environment check failed: {error}\n")

        javac_flags = ["--release", "8"]
        if "java" in args.languages and javac:
            try:
                version = subprocess.run(
                    [str(javac), "-version"], capture_output=True, text=True,
                    timeout=args.timeout,
                )
                if version.returncode:
                    raise RuntimeError((version.stdout + version.stderr).strip())
            except (RuntimeError, OSError, subprocess.TimeoutExpired) as error:
                parser.exit(1, f"Java environment check failed: {error}\n")
            if re.search(r"javac 1\.8\b", version.stdout + version.stderr):
                javac_flags = ["-source", "8", "-target", "8"]

        number = 0
        for number, (guide, line, language, code) in enumerate(examples(args.languages), 1):
            record = {"guide": guide, "line": line, "language": language, "kind": "code block"}
            directory = base / str(number)
            directory.mkdir()
            try:
                if language == "python":
                    if re.search(r"^(?:from rdkit\b|import rdkit\b)", code, re.M) and not rdkit_available:
                        record.update(status="skipped", reason="Optional RDKit dependency is absent")
                    else:
                        source = directory / "example.py"
                        source.write_text(code, encoding="utf-8")
                        record["output"] = run([sys.executable, source], directory, args.timeout)
                elif language == "java":
                    if not java or not javac or not jar or not jar.is_file():
                        raise RuntimeError("Java examples need a JDK and a dependency JAR; use --java-home and --jar")
                    match = re.search(r"public\s+class\s+(\w+)", code)
                    if not match:
                        raise RuntimeError("Java example has no public class")
                    name = match.group(1)
                    source = directory / f"{name}.java"
                    source.write_text(code, encoding="utf-8")
                    run([javac, *javac_flags, "-cp", jar, source], directory, args.timeout)
                    record["output"] = run(
                        [java, "-cp", str(jar) + os.pathsep + str(directory), name],
                        directory, args.timeout,
                    )
                else:
                    if not compiler or not shutil.which(compiler[0]):
                        raise RuntimeError("C++ examples need a C++17 compiler; use --cxx")
                    source = directory / "example.cpp"
                    binary = directory / f"example{suffix}"
                    source.write_text(code, encoding="utf-8")
                    run(
                        [*compiler, "-std=c++17", "-O2", "-I", args.include_dir.resolve(), source, "-o", binary],
                        directory, args.timeout,
                    )
                    record["output"] = run([binary], directory, args.timeout)
                record.setdefault("status", "passed")
            except (RuntimeError, OSError, subprocess.TimeoutExpired) as error:
                record.update(status="failed", reason=str(error))
            report(results, record)

        for path, source in scripts:
            number += 1
            directory = base / str(number)
            directory.mkdir()
            kind = "notebook" if path.suffix == ".ipynb" else "script"
            record = {"guide": path.relative_to(ROOT).as_posix(), "line": 1, "language": "python",
                      "kind": kind}
            try:
                missing = sorted(required_modules(source)) if source is not None else []
                missing = [name for name in missing if name in absent]
                if source is None:
                    record.update(status="skipped", reason="The notebook uses IPython-only syntax")
                elif missing:
                    record.update(status="skipped",
                                  reason=f"Optional {', '.join(missing)} dependency is absent")
                else:
                    if kind == "notebook":
                        target = directory / "notebook.py"
                        target.write_text(source, encoding="utf-8")
                    else:
                        target = path
                    try:
                        record["output"] = run([sys.executable, target], directory, args.timeout)
                    except RuntimeError as error:
                        name = missing_optional(error, absent)
                        if name is None:
                            raise
                        record.update(status="skipped", reason=f"Optional {name} dependency is absent")
                record.setdefault("status", "passed")
            except (RuntimeError, OSError, SyntaxError, subprocess.TimeoutExpired) as error:
                record.update(status="failed", reason=str(error))
            report(results, record)

        for guide, line, code in assistant_loops() if "python" in args.languages else ():
            number += 1
            directory = base / str(number)
            directory.mkdir()
            record = {"guide": guide, "line": line, "language": "python", "kind": "offline tool loop"}
            try:
                source = directory / "loop.py"
                source.write_text(code, encoding="utf-8")
                # Keys the loop reads are set to a placeholder; no request leaves this machine.
                keys = re.findall(r"os\.environ\[\s*[\"'](\w+)[\"']\s*\]", code)
                record["output"] = run([sys.executable, "-c", HARNESS, Path(__file__).resolve(), source],
                                       directory, args.timeout, {key: "offline" for key in keys})
                record["status"] = "passed"
            except (RuntimeError, OSError, subprocess.TimeoutExpired) as error:
                record.update(status="failed", reason=str(error))
            report(results, record)

    counts = {status: sum(item["status"] == status for item in results)
              for status in ("passed", "failed", "skipped")}
    print(f"Examples: {counts['passed']} passed, {counts['failed']} failed, {counts['skipped']} skipped")
    by_language = {
        language: {status: sum(item["language"] == language and item["status"] == status
                              for item in results)
                   for status in counts}
        for language in args.languages
    }
    scope = ("Standalone code blocks, the Python example scripts and notebook, and the assistant "
             "tool loops with fixed offline model replies; CLI, build, installation and "
             "publication commands excluded")
    print(scope)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps({"scope": scope, "environment": environment, "counts": counts,
                        "by_language": by_language, "examples": results}, indent=2) + "\n",
            encoding="utf-8",
        )
    return 1 if counts["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
