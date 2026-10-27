"""Explicit workspace binding independent of where an agent installs the skill."""
import os
from pathlib import Path
from reaction_core import InputError


def bind_workspace(config, project_root=None):
    root = Path(project_root or Path.cwd()).resolve()
    if not root.is_dir():
        raise InputError("selected project root must be an existing directory")
    declared = config.get("project_root")
    if declared is not None and Path(declared).resolve() != root:
        raise InputError("config workspace ownership differs from selected project root")
    if not isinstance(config.get("project_id"), str) or not config["project_id"].strip():
        raise InputError("a nonempty project_id is required for traceability")
    # Do not mutate caller-owned input or couple installation and output locations.
    result = dict(config)
    result["project_root"] = str(root)
    return result


def comsol_install(value=None):
    value = value or os.environ.get("COMSOL_ROOT")
    if not value:
        raise InputError("set COMSOL_ROOT or pass --comsol-root pointing to Multiphysics")
    root = Path(value).resolve()
    for rel in ("bin/win64/comsolbatch.exe", "java/win64/jre/bin/javac.exe"):
        if not (root / rel).is_file():
            raise InputError("COMSOL native executable not found: " + rel)
    return root


def dwsim_install(value=None):
    value = value or os.environ.get("DWSIM_ROOT")
    if not value:
        raise InputError("set DWSIM_ROOT or pass --dwsim-root")
    root = Path(value).resolve()
    if not (root / "DWSIM.Thermodynamics.dll").is_file():
        raise InputError("DWSIM.Thermodynamics.dll not found")
    return root
