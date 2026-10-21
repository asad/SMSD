#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Run the standalone Python, Java and C++ examples in the current guides."""

import argparse
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


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
)
FENCE = re.compile(r"^```(python|java|cpp)[ \t]*\n(.*?)^```[ \t]*$", re.M | re.S)


def run(command, directory, timeout):
    environment = os.environ.copy()
    environment.pop("PYTHONOPTIMIZE", None)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
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


def default_jar():
    metadata = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    version = re.search(r'^version = "([^"]+)"', metadata, re.M).group(1)
    filename = f"smsd-cli-{version}-jar-with-dependencies.jar"
    for directory in (ROOT / "java/cli/target", ROOT / "dist" / f"release-{version}"):
        candidate = directory / filename
        if candidate.is_file():
            return candidate
    return None


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
        if "python" in args.languages:
            try:
                probe = run(
                    [sys.executable, "-c", "import importlib.util, json, sys, smsd; "
                     "print(json.dumps({'python': sys.version.split()[0], "
                     "'smsd': smsd.__version__, 'smsd_path': smsd.__file__, "
                     "'rdkit_available': importlib.util.find_spec('rdkit') is not None}))"],
                    base, args.timeout,
                )
                environment = json.loads(probe)
                rdkit_available = environment["rdkit_available"]
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

        for number, (guide, line, language, code) in enumerate(examples(args.languages), 1):
            record = {"guide": guide, "line": line, "language": language}
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
            results.append(record)
            print(f"{record['status'].upper()}: {guide}:{line} ({language})", flush=True)
            if record.get("reason"):
                print(f"  {record['reason']}", flush=True)

    counts = {status: sum(item["status"] == status for item in results)
              for status in ("passed", "failed", "skipped")}
    print(f"Examples: {counts['passed']} passed, {counts['failed']} failed, {counts['skipped']} skipped")
    by_language = {
        language: {status: sum(item["language"] == language and item["status"] == status
                              for item in results)
                   for status in counts}
        for language in args.languages
    }
    scope = "Standalone code blocks; CLI, build, installation and publication commands excluded"
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
