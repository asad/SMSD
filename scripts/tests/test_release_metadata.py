# SPDX-License-Identifier: Apache-2.0
"""Version updates and the release text share one explicit version."""
import importlib.util
from pathlib import Path
import tempfile

import pytest

SPEC = importlib.util.spec_from_file_location("release_metadata", Path(__file__).parents[1] / "release_metadata.py")
METADATA = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(METADATA)


def test_version_update_preserves_citation_article_date():
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        for name in METADATA.VERSION_FILES:
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('version = "8.0.3"\n')
        (root / "CITATION.cff").write_text('version: 8.0.3\npreferred-citation:\n  date-published: 2026-04-02\n')
        METADATA.set_version(root, "8.1.0", "2026-10-25")
        assert all('8.1.0' in (root / name).read_text() for name in METADATA.VERSION_FILES)
        assert (root / "CITATION.cff").read_text() == (
            'version: 8.1.0\ndate-released: 2026-10-25\npreferred-citation:\n  date-published: 2026-04-02\n')


def test_notes_select_one_release_and_keep_subheadings(tmp_path):
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs/RELEASE_NOTES.md").write_text(
        '# SMSD Pro 8.1.0\n\nCompatible changes.\n\n## Downloads\n\nGitHub.\n\n# SMSD Pro 8.0.3\n\nOlder.\n')
    assert METADATA.release_notes(tmp_path, "8.1.0") == 'Compatible changes.\n\n## Downloads\n\nGitHub.\n'
    with pytest.raises(ValueError, match="missing"):
        METADATA.release_notes(tmp_path, "9.0.0")


@pytest.mark.parametrize("version", ["v8.1.0", "8.1.0\n", "08.1.0", "8.1.0-rc1"])
def test_version_requires_plain_stable_numbers(version):
    with pytest.raises(ValueError):
        METADATA.stable_version(version)
