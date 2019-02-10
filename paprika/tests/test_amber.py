import os
import re
import shutil

import numpy as np
import parmed as pmd

from paprika import amber
from paprika import restraints
from paprika import tleap
from paprika.tests import addons

@pytest.fixture(scope="function", autouse=True)
def clean_files(directory="tmp"):
    # This happens before the test function call
    if os.path.isdir(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    yield
    # This happens after the test function call
    shutil.rmtree(directory)


@addons.using_sander
@addons.using_pmemd_cuda
def test_amber_single_window_gbmin(clean_files):
    # Align PDB to Z-axis
    inputpdb = pmd.load_file('../data/cb6-but/cb6-but-minimized.pdb')
    
    # Distance restraint
    rest1 = restraints.DAT_restraint()
    rest1.continuous_apr = True
    rest1.amber_index = True
    rest1.topology = inputpdb
    rest1.mask1 = ':CB6@O'
    rest1.mask2 = ':BUT@C1'
    rest1.attach['target'] = 4.5
    rest1.attach['fraction_list'] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest1.attach['fc_final'] = 5.0
    rest1.pull['fc'] = rest1.attach['fc_final']
    rest1.pull['target_initial'] = rest1.attach['target']
    rest1.pull['target_final'] = 18.5
    rest1.pull['num_windows'] = 29
    rest1.initialize()
    
    # Angle restraint
    rest2 = restraints.DAT_restraint()
    rest2.continuous_apr = True
    rest2.amber_index = True
    rest2.topology = inputpdb
    rest2.mask1 = ':CB6@O1'
    rest2.mask2 = ':CB6@O'
    rest2.mask3 = ':BUT@C1'
    rest2.attach['target'] = 25.0
    rest2.attach['fraction_list'] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest2.attach['fc_final'] = 50.0
    rest2.pull['fc'] = rest2.attach['fc_final']
    rest2.pull['target_initial'] = rest2.attach['target']
    rest2.pull['target_final'] = rest2.attach['target']
    rest2.pull['num_windows'] = 29
    rest2.initialize()
    
    # Dihedral restraint
    rest3 = restraints.DAT_restraint()
    rest3.continuous_apr = True
    rest3.amber_index = True
    rest3.topology = inputpdb
    rest3.mask1 = ':CB6@O11'
    rest3.mask2 = ':CB6@O1'
    rest3.mask3 = ':CB6@O'
    rest3.mask4 = ':BUT@C1'
    rest3.attach['target'] = -60.0
    rest3.attach['fraction_list'] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest3.attach['fc_final'] = 50.0
    rest3.pull['fc'] = rest3.attach['fc_final']
    rest3.pull['target_initial'] = rest3.attach['target']
    rest3.pull['target_final'] = rest3.attach['target']
    rest3.pull['num_windows'] = 29
    rest3.initialize()

    # Write restraints file for amber
    with open(os.path.join("tmp", 'restraints.in'), 'w') as f:
        for rest in [rest1,rest2,rest3]:
            # Testing just window p001
            f.write(restraints.amber_restraint_line(rest, "p001"))

    # Copy build files for tleap
    files = 'cb6.mol2 cb6.frcmod but.mol2 but.frcmod cb6-but-minimized.pdb'.split()
    for file in files:
        shutil.copy('../data/cb6-but/'+file,path+file)

    # Build prmtop/inpcrd
    sys = tleap.System()
    sys.template_file = './cb6-but/tleap_gb.in'
    sys.output_path = path
    sys.output_prefix = 'vac'
    sys.pbc_type = None
    sys.loadpdb_file = 'cb6-but-minimized.pdb'
    sys.build()

    # Create Simulation
    gbsim = amber.Simulation()
    gbsim.path = path
    gbsim.executable = 'sander'
    gbsim.topology = 'vac.prmtop'
    gbsim.prefix = 'minimize'
    gbsim.inpcrd = 'vac.rst7'
    gbsim.config_gb_min()
    gbsim.cntrl['maxcyc'] = 1
    gbsim.cntrl['ncyc'] = 1 
    gbsim.run()

    # Collect values
    test_vals = []
    with open(path+'minimize.out', 'r') as f:
        filelines = f.readlines()
        for i,line in enumerate(filelines):
            if re.search('^ BOND ', line):
                cols = line.split()
                test_vals.append(float(cols[2])) # BOND
                test_vals.append(float(cols[5])) # ANGLE
                test_vals.append(float(cols[8])) # DIHED
                cols = filelines[i+1].split()
                test_vals.append(float(cols[2])) # VDWAALS
                test_vals.append(float(cols[5])) # EEL
                test_vals.append(float(cols[8])) # EGB
                cols = filelines[i+2].split()
                test_vals.append(float(cols[3])) # 1-4 VDW
                test_vals.append(float(cols[7])) # 1-4 EEL
                test_vals.append(float(cols[10])) # RESTRAINT
                cols = filelines[i+3].split()
                test_vals.append(float(cols[2])) # EAMBER
                cols = filelines[i+4].split()
                test_vals.append(float(cols[4])) # Restraint: Bond
                test_vals.append(float(cols[7])) # Restraint: Angle
                test_vals.append(float(cols[10])) # Restraint: Torsion
                break

    # A bit ugly, but we're here
    test_vals = np.asarray(test_vals)
    # Reference           BOND     ANGLE    DIHED    VDWAALS     EEL        EGB    1-4 VDW    1-4 EEL  RESTRAINT   EAMBER  Bond   Angle  Torsion
    ref_vals = np.array([0.9209, 107.7555, 52.3356, -68.1311, 1331.1371, -127.7541, 5.0130, -2128.7873, 4.5071, -827.5103, 1.214, 2.829, 0.464])

    for i in range(len(ref_vals)):
        assert np.isclose(test_vals[i], ref_vals[i], rtol=0.0, atol=0.0001)

    gbsim.config_gb_md()
    gbsim.prefix = 'md'
    gbsim.inpcrd = 'minimize.rst7'
    gbsim.cntrl['nstlim'] = 1
    gbsim.cntrl['ntpr'] = 1
    gbsim.cntrl['ig'] = 777
    gbsim.run()

    test_vals = []
    with open(path+'md.out', 'r') as f:
        filelines = f.readlines()
        for i,line in enumerate(filelines):
            if re.search('^ NSTEP =        1 ', line):
                cols = line.split()
                test_vals.append(float(cols[8])) # TEMP
                cols = filelines[i+1].split()
                test_vals.append(float(cols[2])) # Etot
                test_vals.append(float(cols[5])) # EKtot
                test_vals.append(float(cols[8])) # EPtot
                cols = filelines[i+2].split()
                test_vals.append(float(cols[2])) # BOND
                test_vals.append(float(cols[5])) # ANGLE
                test_vals.append(float(cols[8])) # DIHED
                cols = filelines[i+3].split()
                test_vals.append(float(cols[3])) # 1-4 NB
                test_vals.append(float(cols[7])) # 1-4 EEL
                test_vals.append(float(cols[10])) # VDWAALS
                cols = filelines[i+4].split()
                test_vals.append(float(cols[2])) # EELEC
                test_vals.append(float(cols[5])) # EGB
                test_vals.append(float(cols[8])) # RESTRAINT
                cols = filelines[i+5].split()
                test_vals.append(float(cols[3])) # EAMBER
                cols = filelines[i+8].split()
                test_vals.append(float(cols[4])) # Restraint: Bond
                test_vals.append(float(cols[7])) # Restraint: Angle
                test_vals.append(float(cols[10])) # Restraint: Torsion
                break

    test_vals = np.asarray(test_vals)
    # Reference           TEMP       Etot    EKtot      EPtot    BOND     ANGLE    DIHED  VDWAALS      EEL        EGB    1-4 VDW     1-4 EEL  RESTRAINT  EAMBER   Bond   Angle  Torsion
    ref_vals = np.array([303.48, -726.5748, 96.4929, -823.0677, 0.8564, 107.7555, 52.3356, 5.0130, -2128.7873, -68.1311, 1331.1371, -127.7541, 4.5071, -827.5748, 1.214, 2.829, 0.464])

    for i in range(len(ref_vals)):
        assert np.isclose(test_vals[i], ref_vals[i], rtol=0.0, atol=0.0001)
