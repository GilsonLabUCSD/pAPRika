"""
Tests the alignment of residues to the z axis.
"""

import unittest
import parmed as pmd
import numpy as np
import logging as log
from paprika.setup.align import *

class TestAlignment(unittest.TestCase):
    
    log.basicConfig(level=log.DEBUG)

    def test_center_mask(self):
        """ Test that the first mask is centered """
        cb6 = pmd.load_file('./cb6-but/vac.pdb')
        aligned_cb6 = align(cb6, ':CB6', ':BUT')
        test_coordinates = check_coordinates(aligned_cb6, ':CB6')
        self.assertTrue(np.allclose(test_coordinates, np.zeros(3)),
                        msg='{}'.format(test_coordinates))

    def test_alignment_after_offset(self):
        """ Test that molecule is properly aligned after random offset. """
        cb6 = pmd.load_file('./cb6-but/vac.pdb')
        random_coordinates = np.random.randint(10) * np.random.rand(1, 3)
        cb6_offset = offset_structure(cb6, random_coordinates)
        aligned_cb6 = align(cb6_offset, ':CB6', ':BUT')
        test_coordinates = check_coordinates(aligned_cb6, ':CB6')
        self.assertTrue(np.allclose(test_coordinates, np.zeros(3)),
                        msg='{}'.format(test_coordinates))


if __name__ == '__main__':
    unittest.main()
    