"""
Tests the solvation of the system using `tleap`.
"""

import unittest
import warnings
import numpy as np
import logging as log
import subprocess as sp
import random as random
import parmed as pmd
import paprika
from paprika.align import *
from paprika.build import *


class TestOpenMM(unittest.TestCase):
    def test_solvation_simple(self):
        """ Test that we can solvate CB6-BUT using default settings. """
        waters = np.random.randint(1000, 10000)
        log.debug('Trying {} waters with default settings...'.format(waters))
        solvate(
            tleapfile='./cb6-but/tleap_solvate.in',
            pdbfile='cb6-but.pdb',
            bufferwater=waters)
        grepped_waters = sp.check_output(
            ["grep -oh 'WAT' ./cb6-but/solvated.prmtop | wc -w"], shell=True)
        self.assertEqual(int(grepped_waters), waters)


if __name__ == '__main__':
    log.debug(f'{paprika.__version__}')
    unittest.main()
