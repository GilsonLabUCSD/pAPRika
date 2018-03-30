"""
Tests tleap tools.
"""

import warnings
import numpy as np
import logging as log
import subprocess as sp
import random as random
import parmed as pmd
import paprika
from paprika.align import *
from paprika.dummy import *
from paprika.tleap import *
import os
import shutil
import sys
import filecmp
import pytest


def rm_solvated_files():
    files = [
        'solvate.pdb', 'solvate.prmtop', 'solvate.rst7', 'tleap_apr_solvate.in', 'leap.log', 'tmp.pdb',
        'cb6-but-dum.tleap.in', 'cb6-but-dum.pdb', 'cb6-but-dum.prmtop', 'cb6-but-dum.rst7', 'tleap.in', 'dm1.mol2', 'dummy.frcmod',
        'solvate.in', 'solvate.tleap.in', 'vac.topo', 'vac.crds', 'full.topo', 'full.crds', 'full.pdb'
    ]
    for f in files:
        if os.path.isfile('./cb6-but/' + f):
            os.remove('./cb6-but/' + f)

@pytest.mark.slow
def test_solvation_simple():
    """ Test that we can solvate CB6-BUT using default settings. """
    waters = np.random.randint(100, 10000)
    #waters = 2922
    log.debug('Trying {} waters with default settings...'.format(waters))
    sys = System()
    sys.template_file = 'cb6-but/tleap_solvate.in'
    sys.output_path = 'cb6-but'
    sys.loadpdb_file = 'cb6-but.pdb'
    sys.buffer_target = waters
    sys.output_prefix = 'solvate'
    sys.build()
    grepped_waters = sp.check_output(["grep -oh 'WAT' ./cb6-but/solvate.prmtop | wc -w"], shell=True)
    assert int(grepped_waters) == waters
    rm_solvated_files()

def test_solvation_octahedron():
    """ Test that we can solvate CB6-BUT with a truncated octahedron. """
    waters = np.random.randint(1000, 10000)
    log.debug('Trying {} waters in a truncated octahedron...'.format(waters))
    sys = System()
    sys.template_file = 'cb6-but/tleap_solvate.in'
    sys.output_path = 'cb6-but'
    sys.loadpdb_file = 'cb6-but.pdb'
    sys.buffer_target = waters
    sys.output_prefix = 'solvate'
    sys.pbc_type = 'octahedral'
    sys.build()
    grepped_waters = sp.check_output(["grep -oh 'WAT' ./cb6-but/solvate.prmtop | wc -w"], shell=True)
    assert int(grepped_waters) == waters
    rm_solvated_files()

def test_solvation_box():
    """ Test that we can solvate CB6-BUT with an isometric box. """
    waters = np.random.randint(1000, 10000)
    log.debug('Trying {} waters in an isometric box...'.format(waters))
    sys = System()
    sys.template_file = 'cb6-but/tleap_solvate.in'
    sys.output_path = 'cb6-but'
    sys.loadpdb_file = 'cb6-but.pdb'
    sys.buffer_target = waters
    sys.output_prefix = 'solvate'
    sys.pbc_type = 'cubic'
    sys.build()
    grepped_waters = sp.check_output(["grep -oh 'WAT' ./cb6-but/solvate.prmtop | wc -w"], shell=True)
    assert int(grepped_waters) == waters
    rm_solvated_files()

@pytest.mark.slow
def test_solvation_spatial_size():
    """ Test that we can solvate CB6-BUT with an buffer size in Angstroms. """
    # this test is a sorta dumb
    random_int = np.random.randint(10, 20)
    random_size = random_int * np.random.random_sample(1) + random_int
    log.debug('Trying buffer size of {} A...'.format(random_size[0]))
    sys = System()
    sys.template_file = 'cb6-but/tleap_solvate.in'
    sys.output_path = 'cb6-but'
    sys.loadpdb_file = 'cb6-but.pdb'
    sys.buffer_target = '{0:1.4f}A'.format(random_size[0])
    sys.output_prefix = 'solvate'
    sys.pbc_type = 'cubic'
    sys.build()
    grepped_waters = sp.check_output(["grep -oh 'WAT' ./cb6-but/solvate.prmtop | wc -w"], shell=True)
    assert int(grepped_waters) == sys.target_waters
    rm_solvated_files()

@pytest.mark.slow
def test_solvation_potassium_control():
    """ Test there is no potassium by default. A negative control. """
    waters = np.random.randint(1000, 10000)
    log.debug('Trying {} waters with potassium...'.format(waters))
    sys = System()
    sys.template_file = 'cb6-but/tleap_solvate.in'
    sys.output_path = 'cb6-but'
    sys.loadpdb_file = 'cb6-but.pdb'
    sys.buffer_target = waters
    sys.output_prefix = 'solvate'
    sys.counter_cation = 'K+'
    sys.build()
    potassium = sp.check_output(["grep -oh 'K+' ./cb6-but/solvate.prmtop | wc -w"], shell=True)
    assert int(potassium) == 0
    rm_solvated_files()

@pytest.mark.slow
def test_solvation_with_additional_ions():
    """ Test that we can solvate CB6-BUT with additional ions. """
    waters = np.random.randint(1000, 10000)
    cations = ['LI', 'Na+', 'K+', 'RB', 'CS']
    anions = ['F', 'Cl-', 'BR', 'IOD']
    n_cations = np.random.randint(1, 10)
    n_anions = np.random.randint(1, 10)
    random_cation = random.choice(cations)
    random_anion = random.choice(anions)
    log.debug('Trying {} waters with additional ions...'.format(waters))
    sys = System()
    sys.template_file = 'cb6-but/tleap_solvate.in'
    sys.output_path = 'cb6-but'
    sys.loadpdb_file = 'cb6-but.pdb'
    sys.buffer_target = waters
    sys.output_prefix = 'solvate'
    sys.neutralize = False
    sys.add_ions = [random_cation, n_cations, random_anion, n_anions]
    sys.build()
    # These should come in the RESIDUE_LABEL region of the prmtop and be before all the water.
    cation_number = sp.check_output(
        ["grep -A 99 RESIDUE_LABEL ./cb6-but/solvate.prmtop | " + "grep -oh '{} ' | wc -w".format(random_cation)],
        shell=True)
    anion_number = sp.check_output(
        ["grep -A 99 RESIDUE_LABEL ./cb6-but/solvate.prmtop | " + "grep -oh '{} ' | wc -w".format(random_anion)],
        shell=True)
    # Have to think about what to do here...
    log.debug('Expecting...')
    log.debug('cation = {}\tn_cations={}'.format(random_cation, n_cations))
    log.debug('anion  = {}\t n_anions={}'.format(random_anion, n_anions))
    log.debug('Found...')
    log.debug('             n_cations={}'.format(cation_number))
    log.debug('              n_anions={}'.format(anion_number))

    assert int(cation_number) == n_cations and int(anion_number) == n_anions
    rm_solvated_files()

def test_solvation_by_M_and_m():
    """ Test that we can solvate CB6-BUT through molarity and molality. """
    log.debug('Trying 10 A buffer with 150 mM NaCl...')
    sys = System()
    sys.template_file = 'cb6-but/tleap_solvate.in'
    sys.output_path = 'cb6-but'
    sys.loadpdb_file = 'cb6-but.pdb'
    sys.buffer_target = '10A'
    sys.output_prefix = 'solvate'
    sys.neutralize = False
    sys.pbc_type = 'rectangular'
    sys.add_ions = ['NA', '0.150M', 'CL', '0.150M', 'K', '0.100m', 'BR', '0.100m']
    sys.build()

    # Molarity Check
    obs_num_na = sp.check_output(
        ["grep -A 99 RESIDUE_LABEL ./cb6-but/solvate.prmtop | " + "grep -oh 'NA ' | wc -w"], shell=True)
    obs_num_cl = sp.check_output(
        ["grep -A 99 RESIDUE_LABEL ./cb6-but/solvate.prmtop | " + "grep -oh 'CL ' | wc -w"], shell=True)

    volume = sys.get_volume()
    volume_in_liters = volume * ANGSTROM_CUBED_TO_LITERS
    calc_num_na = np.ceil((6.022 * 10**23) * (0.150) * volume_in_liters)
    calc_num_cl = np.ceil((6.022 * 10**23) * (0.150) * volume_in_liters)
    assert int(obs_num_na) == calc_num_na and int(obs_num_cl) == calc_num_cl

    # Molality Check
    obs_num_k = sp.check_output(
        ["grep -A 99 RESIDUE_LABEL ./cb6-but/solvate.prmtop | " + "grep -oh 'K ' | wc -w"], shell=True)
    obs_num_br = sp.check_output(
        ["grep -A 99 RESIDUE_LABEL ./cb6-but/solvate.prmtop | " + "grep -oh 'BR ' | wc -w"], shell=True)
    calc_num_waters = sys.count_residues()['WAT']
    calc_num_k = np.ceil(0.100 * calc_num_waters * 0.018)
    calc_num_br = np.ceil(0.100 * calc_num_waters * 0.018)
    assert int(obs_num_k) == calc_num_k and int(obs_num_br) == calc_num_br
    rm_solvated_files()

@pytest.mark.slow
def test_alignment_workflow():
    """ Test that we can solvate CB6-BUT after alignment. """
    cb6 = pmd.load_file('./cb6-but/cb6-but-notcentered.pdb')
    zalign(cb6, ':CB6', ':BUT', save=True, filename='./cb6-but/tmp.pdb')
    waters = np.random.randint(1000, 10000)
    sys = System()
    sys.template_file = 'cb6-but/tleap_solvate.in'
    sys.output_path = 'cb6-but'
    sys.loadpdb_file = 'tmp.pdb'
    sys.buffer_target = waters
    sys.output_prefix = 'solvate'
    sys.build()
    log.debug('Trying {} waters after alignment...'.format(waters))
    grepped_waters = sp.check_output(["grep -oh 'WAT' ./cb6-but/solvate.prmtop | wc -w"], shell=True)
    assert int(grepped_waters) == waters
    rm_solvated_files()

def test_add_dummy():
    """ Test that dummy atoms get added correctly """
    hostguest = pmd.load_file('cb6-but/cb6-but-notcentered.pdb', structure=True)
    hostguest = zalign(hostguest, ':BUT@C', ':BUT@C3', save=False)
    hostguest = add_dummy(hostguest, residue_name='DM1', z=-11.000, y=2.000, x=-1.500)
    hostguest.write_pdb('cb6-but/cb6-but-dum.pdb', renumber=False)
    with open('cb6-but/cb6-but-dum.pdb', 'r') as f:
        lines = f.readlines()
        test_line1 = lines[123].rstrip()
        test_line2 = lines[124].rstrip()
    ref_line1 = 'TER     123      BUT     2'
    ref_line2 = 'HETATM  123 DUM  DM1     3      -1.500   2.000 -11.000  0.00  0.00          PB'
    assert ref_line1 == test_line1
    assert ref_line2 == test_line2

    write_dummy_frcmod(path='./cb6-but/')
    write_dummy_mol2(path='./cb6-but/', filename='dm1.mol2', residue_name='DM1')
    sys = System()
    sys.template_lines = [
        'source leaprc.gaff',
        'loadamberparams cb6.frcmod',
        'CB6 = loadmol2 cb6.mol2',
        'loadamberparams but.frcmod',
        'BUT = loadmol2 but.mol2',
        'loadamberparams dummy.frcmod',
        'DM1 = loadmol2 dm1.mol2',
        'model = loadpdb cb6-but-dum.pdb'
    ]
    sys.output_path = 'cb6-but'
    sys.output_prefix = 'cb6-but-dum'
    sys.pbc_type = None
    sys.build()

    assert filecmp.cmp('cb6-but/REF_cb6-but-dum.rst7', 'cb6-but/cb6-but-dum.rst7', shallow=False)
    rm_solvated_files()


