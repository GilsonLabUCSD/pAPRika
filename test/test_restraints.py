"""
Tests the restraints utilities.
"""

import unittest
import warnings
import numpy as np
import logging as log
import subprocess as sp
import parmed as pmd
import pytest
import paprika
from paprika.restraints import *

def test_DAT_restraint():
    rest1 = restraints.DAT_restraint()
    rest1.continuous_apr = True
    rest1.structure_file = '../cb6-but/cb6-but-notcentered.pdb'
    rest1.mask1 = ':CB6@O,O2,O4,O6,O8,O10'
    rest1.mask2 = ':BUT@C*'
    rest1.attach['target'] = 3.0
    rest1.attach['fraction_list'] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest1.attach['fc_final'] = 5.0
    rest1.pull['fc'] = rest1.attach['fc_final']
    rest1.pull['target_initial'] = rest1.attach['target']
    rest1.pull['target_final'] = 10.0
    rest1.pull['num_windows'] = 11
    rest1.initialize()

    assert rest1.index1 == [13, 31, 49, 67, 85, 103]
    assert rest1.index2 == [109, 113, 115, 119]
    assert rest1.phase['attach']['force_constants'] == [0.0, 0.2, 0.905, 2.48, 5.0]
    assert rest1.phase['attach']['targets'] == [3.0, 3.0, 3.0, 3.0, 3.0]
    assert rest1.phase['pull']['force_constants'] == [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    assert rest1.phase['pull']['targets'] == [  3.    3.7   4.4   5.1   5.8   6.5   7.2   7.9   8.6   9.3  10. ]
    assert rest1.phase['release']['force_constants'] == [ 5.  5.  5.  5.  5.]
    assert rest1.phase['release']['targets'] == [3.0, 3.0, 3.0, 3.0, 3.0]
