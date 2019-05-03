"""
Decorators and wrappers for testing.
"""

import pkgutil
import shutil

import pytest


def _plugin_import(plug):
    plug_spec = pkgutil.find_loader(plug)
    if plug_spec is None:
        return False
    else:
        return True


def _find_executable(exe):
    if shutil.which(exe) is not None:
        return True
    else:
        return False


_import_message = (
    "Cannot import module {}. Install package if necessary and add to PYTHONPATH"
)
_exe_message = "Cannot detect executable {}. Install package if necessary."

using_openmm = pytest.mark.skipif(
    _plugin_import("simtk.openmm") is False, reason=_import_message.format("OpenMM")
)
using_sander = pytest.mark.skipif(
    _find_executable("sander") is False, reason=_import_message.format("sander")
)
using_pmemd_cuda = pytest.mark.skipif(
    _find_executable("pmemd.cuda") is False, reason=_import_message.format("pmemd.cuda")
)
using_tleap = pytest.mark.skipif(
    _find_executable("tleap") is False, reason=_import_message.format("tLEaP")
)

using_parmed = pytest.mark.skipif(
    _plugin_import("parmed") is False, reason=_import_message.format("ParmEd")
)
