"""
Unit and regression test for the paprika package.
"""

import sys

# Import package, test suite, and other packages as needed
import paprika
import pytest


def test_paprika_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "paprika" in sys.modules
