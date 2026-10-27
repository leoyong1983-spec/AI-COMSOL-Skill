"""Regression checks for relocation, ownership and artifact containment."""
import json
import unittest
from pathlib import Path
from unittest.mock import patch
from reaction_core import InputError, owned_path
from workspace import bind_workspace, comsol_install
from electrosynthesis import load

ASSETS = Path(__file__).resolve().parents[1] / "assets"


class WorkspaceTests(unittest.TestCase):
    def test_input_does_not_depend_on_skill_installation_location(self):
        c = json.loads((ASSETS / "literature_eg_equivalent.json").read_text(encoding="utf-8"))
        model = load(ASSETS / "literature_eg_equivalent.json", Path.cwd())
        self.assertEqual(Path(model.c["project_root"]), Path.cwd().resolve())
        self.assertNotIn("project_root", c)

    def test_declared_different_project_is_blocked(self):
        with self.assertRaisesRegex(InputError, "ownership"):
            bind_workspace({"project_id": "demo", "project_root": str(Path.cwd().parent)}, Path.cwd())

    def test_relative_artifacts_follow_selected_workspace(self):
        root = Path.cwd().resolve()
        self.assertEqual(owned_path(root, "AGENTCODING/07_handoff/new", "07_handoff"), root / "AGENTCODING/07_handoff/new")
        with self.assertRaises(InputError):
            owned_path(root, "AGENTCODING/07_handoff/../../escape", "07_handoff")

    def test_wrong_artifact_stage_is_blocked(self):
        with self.assertRaises(InputError):
            owned_path(Path.cwd(), "AGENTCODING/04_work/model.mph", "07_handoff")

    def test_missing_installation_has_actionable_error(self):
        with patch.dict("os.environ", {}, clear=True):
            with self.assertRaisesRegex(InputError, "COMSOL_ROOT"):
                comsol_install()


if __name__ == "__main__":
    unittest.main(verbosity=2)
