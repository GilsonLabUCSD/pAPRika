"""
Decorators and wrappers for testing.
"""

import pkgutil
import pytest
import subprocess as sp

def _plugin_import(plug):
    plug_spec = pkgutil.find_loader(plug)
    if plug_spec is None:
        return False
    else:
        return True

def _find_executible(exe):
    # Switch to shutil.which in Py3.3
    try:
        sp.check_out(["which", exe])
        return True
    except sp.CalledProcessError:
        return False


_import_message = "Not detecting module {}. Install package if necessary and add to envvar PYTHONPATH"
_exe_message = "Not detecting executable {}. Install package if necessary."

using_openmm = pytest.mark.skipif(_plugin_import("simtk.openmm") is False, reason=_import_message.format('OpenMM'))
using_amber = pytest.mark.skipif(_find_executible("amber") is False, reason=_import_message.format('Amber'))
using_tleap = pytest.mark.skipif(_find_executible("tleap") is False, reason=_import_message.format('tLeap'))

