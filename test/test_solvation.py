"""
Tests the solvation of the system using `tleap`.
"""

import unittest
import parmed as pmd
import numpy as np
import logging as log
from paprika.setup.align import *
from paprika.setup.solvate import *

class TestSolvate(unittest.TestCase):
    
    logger = log.getLogger()
    logger.setLevel(log.DEBUG)

    def test_solvation_simple(self):
        """ Test that we can solvate CB6-BUT using default settings. """
        waters = np.random.randint(1000, 10000)
        # log.debug('Trying to solvate with {} waters...'.format(waters))
        solvate(tleapfile='./cb6-but/tleap.in', pdbfile='cb6-but.pdb',
                bufferwater=waters)
        self.assertEqual(countresidues(filename='tleap_apr_solvate.in', directory='./cb6-but/', returnlist='WAT'), waters, msg='{}'.format(waters))

     def test_alignment_workflow(self):
         """ Test that we can solvate CB6-BUT after alignment. """
         cb6 = pmd.load_file('./cb6-but/vac.pdb')
         align(cb6, ':CB6', ':BUT', save=True, filename='./cb6-but/tmp.pdb')
         waters = np.random.randint(1000, 10000)
         print('Trying to solvate with {} waters...'.format(waters))
         solvate(tleapfile='./cb6-but/tleap.in', pdbfile='tmp.pdb',
                 bufferwater=waters)
         self.assertEqual(countresidues(filename='tleap_apr_solvate.in', directory='./cb6-but/',returnlist='WAT'), waters, msg='{}'.format(waters))

if __name__ == '__main__':
    unittest.main()
    
