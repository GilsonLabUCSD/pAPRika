"""
Unit and regression test for the paprika package.
"""

# Import package, test suite, and other packages as needed
import paprika
import pytest
import sys


def test_paprika_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "paprika" in sys.modules
