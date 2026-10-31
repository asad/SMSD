# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Offline checks of the example checker; no built smsd, model service or network is needed."""
import importlib.util
import json
from pathlib import Path
import re
import sys
import tempfile
import textwrap
import unittest

SCRIPTS = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CHECK = load("check_examples_under_test", SCRIPTS / "check_examples.py")

# A stand-in smsd.tool_definitions with the published call shapes and argument checks.
FAKE_TOOLS = '''
_SCHEMA = {"type": "object", "additionalProperties": False, "required": ["query_smiles", "target_smiles"],
           "properties": {"query_smiles": {"type": "string"}, "target_smiles": {"type": "string"}}}


def get_tool_definitions(provider="mcp"):
    base = {"name": "find_mcs", "description": "Find one common graph."}
    if provider == "input-schema":
        return [dict(base, input_schema=_SCHEMA)]
    if provider == "function":
        return [dict(base, type="function", parameters=_SCHEMA, strict=False)]
    if provider == "chat-function":
        return [{"type": "function", "function": dict(base, parameters=_SCHEMA, strict=False)}]
    return [dict(base, inputSchema=_SCHEMA)]


def call_tool(name, arguments):
    if name != "find_mcs":
        raise ValueError("Unknown tool")
    if any(key not in arguments for key in _SCHEMA["required"]):
        raise ValueError("Missing required tool argument")
    return {"found": True, "status": "COMPLETE", "optimal": True, "complete": True, "score": 6,
            "matches": [{"atom_count": 6}]}
'''


class ImportRules(unittest.TestCase):
    def test_fallback_imports_are_optional_and_exiting_imports_are_required(self):
        fallback = "try:\n    import smsd_cuda\nexcept ImportError:\n    smsd_cuda = None\n"
        self.assertEqual(CHECK.required_modules(fallback), set())
        exiting = ("import sys\ntry:\n    from mcp import ClientSession\n"
                   "except ImportError:\n    sys.exit('needs the MCP SDK')\n")
        self.assertEqual(CHECK.required_modules(exiting), {"mcp"})
        raising = "try:\n    import numpy\nexcept ImportError:\n    raise\n"
        self.assertEqual(CHECK.required_modules(raising), {"numpy"})

    def test_standard_library_smsd_and_conditional_imports_are_not_required(self):
        source = ("import json\nimport smsd\nfrom smsd.corpus import FingerprintCorpus\n"
                  "if False:\n    from google.colab import files\nfrom rdkit import Chem\n")
        self.assertEqual(CHECK.required_modules(source), {"rdkit"})

    def test_only_a_final_import_error_naming_an_absent_module_skips(self):
        error = RuntimeError("Command exited 1: Traceback\n  File x\n"
                             "ImportError: FingerprintCorpus requires NumPy; install numpy>=1.23")
        self.assertEqual(CHECK.missing_optional(error, {"numpy"}), "numpy")
        self.assertIsNone(CHECK.missing_optional(error, set()))
        native = RuntimeError("Command exited 1:\nModuleNotFoundError: No module named 'smsd._smsd'")
        self.assertIsNone(CHECK.missing_optional(native, {"numpy"}))
        later = RuntimeError("Command exited 1:\nImportError: numpy\nAssertionError")
        self.assertIsNone(CHECK.missing_optional(later, {"numpy"}))


class ExampleDiscovery(unittest.TestCase):
    def test_every_example_script_and_notebook_is_run(self):
        found = {path.name: source for path, source in CHECK.example_scripts()}
        directory = SCRIPTS.parent / "python/examples"
        expected = {path.name for path in directory.iterdir() if path.suffix in (".py", ".ipynb")}
        self.assertEqual(set(found), expected)
        for name, source in found.items():
            self.assertIsNotNone(source, name)
            compile(source, name, "exec")

    def test_optional_example_scripts_are_classified_when_present(self):
        directory = SCRIPTS.parent / "python/examples"
        expectations = {"gpu_screening.py": set(), "agent_workflow.py": set(),
                        "mcp_stdio_client.py": {"mcp"}, "feature_adapters.py": set()}
        for name, required in expectations.items():
            path = directory / name
            if path.is_file():
                self.assertEqual(CHECK.required_modules(path.read_text(encoding="utf-8")), required, name)

    def test_notebook_cells_are_joined_and_ipython_syntax_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "example.ipynb"
            cells = [{"cell_type": "markdown", "source": ["# Title"]},
                     {"cell_type": "code", "source": ["value = 1\n", "print(value)"]}]
            path.write_text(json.dumps({"cells": cells}))
            self.assertIn("print(value)", CHECK.notebook_source(path))
            cells.append({"cell_type": "code", "source": ["%timeit value"]})
            path.write_text(json.dumps({"cells": cells}))
            self.assertIsNone(CHECK.notebook_source(path))


class OfflineToolLoops(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        package = self.root / "fake/smsd"
        package.mkdir(parents=True)
        (package / "__init__.py").write_text('__version__ = "0"\n')
        (package / "tool_definitions.py").write_text(textwrap.dedent(FAKE_TOOLS))

    def tearDown(self):
        self.directory.cleanup()

    def drive(self, code, name="loop.py"):
        source = self.root / name
        source.write_text(code, encoding="utf-8")
        keys = {key: "offline" for key in re.findall(r"os\.environ\[\s*[\"'](\w+)[\"']\s*\]", code)}
        keys["PYTHONPATH"] = str(self.root / "fake")
        return CHECK.run([sys.executable, "-c", CHECK.HARNESS, SCRIPTS / "check_examples.py", source],
                         self.root, 120, keys)

    def test_each_published_loop_returns_results_and_errors_to_the_model(self):
        loops = list(CHECK.assistant_loops())
        self.assertTrue(loops)
        styles = set()
        for guide, line, code in loops:
            with self.subTest(line=line):
                summary = json.loads(self.drive(code).splitlines()[-1])
                self.assertEqual(summary["tool_results"], 2)
                styles.add(summary["reply_style"])
        self.assertEqual(styles, set(CHECK.ENDPOINTS.values()))

    def test_a_loop_that_lets_tool_errors_escape_fails(self):
        _, _, code = next(iter(CHECK.assistant_loops()))
        broken = code.replace("except (ValueError, ImportError) as error:", "except ImportError as error:")
        self.assertNotEqual(broken, code)
        with self.assertRaisesRegex(RuntimeError, "Missing required tool argument"):
            self.drive(broken)

    def test_a_loop_that_never_asks_the_model_fails(self):
        code = "import json\nfrom smsd.tool_definitions import call_tool\nprint(json.dumps({}))\n"
        with self.assertRaisesRegex(RuntimeError, "no model client"):
            self.drive(code)


if __name__ == "__main__":
    unittest.main()
