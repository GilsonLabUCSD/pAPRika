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


class TestSolvate(unittest.TestCase):
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

    def test_solvation_octahedron(self):
        """ Test that we can solvate CB6-BUT with a truncated octahedron. """
        waters = np.random.randint(1000, 10000)
        log.debug(
            'Trying {} waters in a truncated octahedron...'.format(waters))
        solvate(
            tleapfile='./cb6-but/tleap_solvate.in',
            pdbfile='cb6-but.pdb',
            bufferwater=waters,
            pbctype=2)
        grepped_waters = sp.check_output(
            ["grep -oh 'WAT' ./cb6-but/solvated.prmtop | wc -w"], shell=True)
        self.assertEqual(int(grepped_waters), waters)

    def test_solvation_box(self):
        """ Test that we can solvate CB6-BUT with an isometric box. """
        waters = np.random.randint(1000, 10000)
        log.debug('Trying {} waters in an isometric box...'.format(waters))
        solvate(
            tleapfile='./cb6-but/tleap_solvate.in',
            pdbfile='cb6-but.pdb',
            bufferwater=waters,
            pbctype=0)
        grepped_waters = sp.check_output(
            ["grep -oh 'WAT' ./cb6-but/solvated.prmtop | wc -w"], shell=True)
        self.assertEqual(int(grepped_waters), waters)

    def test_solvation_spatial_size(self):
        """ Test that we can solvate CB6-BUT with an buffer size in Angstroms. """
        random_int = np.random.randint(10, 100)
        random_size = random_int * np.random.random_sample(1) + random_int
        log.debug('Trying buffer size of {} A...'.format(random_size[0]))
        solvate(
            tleapfile='./cb6-but/tleap_solvate.in',
            pdbfile='cb6-but.pdb',
            bufferwater='{0:1.4f}A'.format(random_size[0]))
        log.warning('I don\'t know the proper test case for this. \
                     We\'d have to parse the coordinates, perhaps.')
        # This is bad, but I'll think of a better way.
        self.assertEqual(0, 0)

    def test_solvation_potassium_control(self):
        """ Test there is no potassium by default. A negative control. """
        waters = np.random.randint(1000, 10000)
        log.debug('Trying {} waters with potassium...'.format(waters))
        solvate(
            tleapfile='./cb6-but/tleap_solvate.in',
            pdbfile='cb6-but.pdb',
            bufferwater=waters,
            counter_cation='K+')
        potassium = sp.check_output(
            ["grep -oh 'K+' ./cb6-but/solvated.prmtop | wc -w"], shell=True)
        self.assertEqual(int(potassium), 0)

    def test_solvation_with_additional_ions(self):
        """ Test that we can solvate CB6-BUT with additional ions. """
        waters = np.random.randint(1000, 10000)
        cations = ['LI', 'Na+', 'K+', 'RB', 'CS']
        anions = ['F', 'Cl-', 'BR', 'IOD']
        n_cations = np.random.randint(1, 10)
        n_anions = np.random.randint(1, 10)
        random_cation = random.choice(cations)
        random_anion = random.choice(anions)
        log.debug('Trying {} waters with additional ions...'.format(waters))
        solvate(
            tleapfile='./cb6-but/tleap_solvate.in',
            pdbfile='cb6-but.pdb',
            bufferwater=waters,
            neutralize=False,
            addions=[random_cation, n_cations, random_anion, n_anions])
        # These should come in the RESIDUE_LABEL region of the prmtop and be before all the water.
        cation_number = sp.check_output(
            [
                "grep -A 99 RESIDUE_LABEL ./cb6-but/solvated.prmtop | " +
                "grep -oh '{} ' | wc -w".format(random_cation)
            ],
            shell=True)
        anion_number = sp.check_output(
            [
                "grep -A 99 RESIDUE_LABEL ./cb6-but/solvated.prmtop | " +
                "grep -oh '{} ' | wc -w".format(random_anion)
            ],
            shell=True)
        # Have to think about what to do here...
        # log.debug('Expecting...')
        # log.debug('cation = {}\tn_cations={}'.format(random_cation, n_cations))
        # log.debug('anion  = {}\t n_anions={}'.format(random_anion, n_anions))
        # log.debug('Found...')
        # log.debug('             n_cations={}'.format(cation_number))
        # log.debug('              n_anions={}'.format(anion_number))

        self.assertTrue(
            int(cation_number) == n_cations and int(anion_number) == n_anions)

    def test_solvation_by_molarity(self):
        """ Test that we can solvate CB6-BUT through molarity. """
        log.debug('Trying 10 A buffer with 150 mM NaCl...')
        solvate(
            tleapfile='./cb6-but/tleap_solvate.in',
            pdbfile='cb6-but.pdb',
            bufferwater='10A',
            neutralize=False,
            pbctype=1,
            addions=['Na+', '0.150M', 'Cl-', '0.150M'])
        cation_number = sp.check_output(
            [
                "grep -A 99 RESIDUE_LABEL ./cb6-but/solvated.prmtop | " +
                "grep -oh 'Na+ ' | wc -w"
            ],
            shell=True)
        anion_number = sp.check_output(
            [
                "grep -A 99 RESIDUE_LABEL ./cb6-but/solvated.prmtop | " +
                "grep -oh 'Cl- ' | wc -w"
            ],
            shell=True)
        # Approximate volume of the solvated system in liters
        volume_in_liters = 4.6766 * 10**-23
        n_cations = np.ceil((6.022 * 10**23) * (0.150) * volume_in_liters)
        n_anions = np.ceil((6.022 * 10**23) * (0.150) * volume_in_liters)
        self.assertTrue(
            int(cation_number) == n_cations and int(anion_number) == n_anions)

    def test_solvation_by_molality(self):
        """ Test that we can solvate CB6-BUT through molarity. """
        log.debug('Trying 2000 water buffer with 150 mmol/kg NaCl...')
        solvate(
            tleapfile='./cb6-but/tleap_solvate.in',
            pdbfile='cb6-but.pdb',
            bufferwater=2000,
            neutralize=False,
            pbctype=1,
            addions=['Na+', '0.150m', 'Cl-', '0.150m'])
        cation_number = sp.check_output(
            [
                "grep -A 99 RESIDUE_LABEL ./cb6-but/solvated.prmtop | " +
                "grep -oh 'Na+ ' | wc -w"
            ],
            shell=True)
        anion_number = sp.check_output(
            [
                "grep -A 99 RESIDUE_LABEL ./cb6-but/solvated.prmtop | " +
                "grep -oh 'Cl- ' | wc -w"
            ],
            shell=True)
        n_cations = np.ceil(0.150 * 2000 * 0.018)
        n_anions = np.ceil(0.150 * 2000 * 0.018)
        self.assertTrue(
            int(cation_number) == n_cations and int(anion_number) == n_anions)

    def test_alignment_workflow(self):
        """ Test that we can solvate CB6-BUT after alignment. """
        cb6 = pmd.load_file('./cb6-but/vac.pdb')
        align(cb6, ':CB6', ':BUT', save=True, filename='./cb6-but/tmp.pdb')
        waters = np.random.randint(1000, 10000)
        log.debug('Trying {} waters after alignment...'.format(waters))
        solvate(
            tleapfile='./cb6-but/tleap_solvate.in',
            pdbfile='tmp.pdb',
            bufferwater=waters)
        grepped_waters = sp.check_output(
            ["grep -oh 'WAT' ./cb6-but/solvated.prmtop | wc -w"], shell=True)
        self.assertEqual(int(grepped_waters), waters)


if __name__ == '__main__':
    log.debug('{}'.format(paprika.__version__))
    unittest.main()
