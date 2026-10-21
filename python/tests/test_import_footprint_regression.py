# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Importing the package, and using its core, loads no optional toolkit.

RDKit, NumPy and Open Babel are imported only by the functions that need them,
so ``import smsd`` stays light and the core works where they are not installed.
Each check runs in a fresh interpreter, because this test session may already
have loaded them.
"""
import importlib.util
import json
import os
import subprocess
import sys
import textwrap

import pytest
import smsd

OPTIONAL_TOOLKITS = {"rdkit", "numpy", "openbabel", "pybel", "pandas", "PIL"}


def _loaded_by(source, tmp_path):
    """The top-level modules that ``source`` adds to a fresh interpreter's ``sys.modules``."""
    script = textwrap.dedent("""
        import json, sys
        before = set(sys.modules)
    """) + textwrap.dedent(source) + textwrap.dedent("""
        print(json.dumps(sorted({name.split(".")[0] for name in set(sys.modules) - before})))
    """)
    # Run against the package this session imported, whether installed or in a source tree.
    environment = dict(os.environ)
    package_root = os.path.dirname(os.path.dirname(os.path.abspath(smsd.__file__)))
    environment["PYTHONPATH"] = os.pathsep.join(filter(None, [package_root, environment.get("PYTHONPATH")]))
    result = subprocess.run([sys.executable, "-c", script], cwd=tmp_path, env=environment,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    return set(json.loads(result.stdout.splitlines()[-1]))


def test_import_loads_no_optional_toolkit(tmp_path):
    loaded = _loaded_by("import smsd", tmp_path)
    assert "smsd" in loaded
    assert not loaded & OPTIONAL_TOOLKITS


def test_import_leaves_pathlib_to_the_file_functions(tmp_path):
    """``pathlib``, and the ``glob`` and ``fnmatch`` it brings, are imported only to read or write a file."""
    loaded = _loaded_by("""
        import sys
        preloaded = "pathlib" in sys.modules
        import smsd
        assert preloaded or "pathlib" not in sys.modules, "import smsd loaded pathlib"
    """, tmp_path)
    assert "smsd" in loaded
    assert not loaded & {"pathlib", "glob", "fnmatch"}


def test_core_workflow_loads_no_optional_toolkit(tmp_path):
    loaded = _loaded_by("""
        import smsd
        phenol = smsd.parse_smiles("c1ccc(O)cc1")
        assert smsd.Substructure.find("c1ccccc1").match(phenol).found
        assert smsd.MCS.find(phenol).match("c1ccc(N)cc1").found
        assert smsd.Substructure.smarts("[OX2H]").match(phenol).found
        reordered = smsd.parse_smiles("Oc1ccccc1")
        assert smsd.tanimoto_coefficient(smsd.circular_fingerprint(phenol), smsd.circular_fingerprint(reordered)) == 1.0
        assert smsd.canonical_smiles("Oc1ccccc1") == smsd.canonical_smiles(phenol)
        assert smsd.generate_coords_2d(phenol)
        assert smsd.depict_svg(phenol).lstrip().startswith("<")
        smsd.write_molfile(phenol, "phenol.mol")
        assert smsd.read_mol_file("phenol.mol").n == phenol.n
    """, tmp_path)
    assert not loaded & OPTIONAL_TOOLKITS


@pytest.mark.skipif(importlib.util.find_spec("rdkit") is None, reason="RDKit is not installed")
def test_rdkit_is_imported_when_a_conversion_needs_it(tmp_path):
    loaded = _loaded_by("""
        import smsd
        assert smsd.from_rdkit(smsd.to_rdkit("CCO")).n == 3
    """, tmp_path)
    assert "rdkit" in loaded
