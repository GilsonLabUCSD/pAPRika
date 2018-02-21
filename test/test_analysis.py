import parmed as pmd
import pytraj as pt
import numpy as np
import paprika.restraints
import paprika.analysis

def test_fe_calc():

    inputpdb = pmd.load_file('cb6_but_gb_ref_data/vac.pdb')
    
    # Distance restraint
    rest1 = paprika.restraints.DAT_restraint()
    rest1.continuous_apr = True
    rest1.amber = True
    rest1.topology = inputpdb
    rest1.mask1 = ':CB6@O'
    rest1.mask2 = ':BUT@C1'
    rest1.attach['target'] = 4.5
    rest1.attach['fraction_list'] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest1.attach['fc_final'] = 5.0
    rest1.pull['fc'] = rest1.attach['fc_final']
    rest1.pull['target_initial'] = rest1.attach['target']
    rest1.pull['target_final'] = 18.5
    rest1.pull['num_windows'] = 19
    rest1.initialize()
    
    # Angle restraint
    rest2 = paprika.restraints.DAT_restraint()
    rest2.continuous_apr = True
    rest2.amber = True
    rest2.topology = inputpdb
    rest2.mask1 = ':CB6@O1'
    rest2.mask2 = ':CB6@O'
    rest2.mask3 = ':BUT@C1'
    rest2.attach['target'] = 8.0
    rest2.attach['fraction_list'] = [0.00, 0.04, 0.181, 0.496, 1.000]
    rest2.attach['fc_final'] = 50.0
    rest2.pull['fc'] = rest2.attach['fc_final']
    rest2.pull['target_initial'] = rest2.attach['target']
    rest2.pull['target_final'] = rest2.attach['target']
    rest2.pull['num_windows'] = 19
    rest2.initialize()
    
    # Dihedral restraint
    rest3 = paprika.restraints.DAT_restraint()
    rest3.continuous_apr = True
    rest3.amber = True
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
    rest3.pull['num_windows'] = 19
    rest3.initialize()
    
    # Create window directories
    window_list = paprika.restraints.create_window_list([rest1,rest2,rest3])

    # Phase abbreviations
    phase_dict = {'a': 'attach', 'p': 'pull', 'r': 'release'}
    
    # Restraint data
    rest_dat = {'attach': [], 'pull': [], 'release': []}
    
    # Loop through windows
    for i,window in enumerate(window_list):
        phase = phase_dict[window[0]]
        rest_dat[phase].append([])
        traj = pt.load('cb6_but_gb_ref_data/'+window+'/md.nc', 'cb6_but_gb_ref_data/vac.prmtop')
        for rest in [rest1,rest2,rest3]:
            if rest.mask1 and rest.mask2 and not rest.mask3 and not rest.mask4:
                dist = pt.distance(traj, ' '.join([rest.mask1, rest.mask2]))
                rest_dat[phase][-1].append(dist)
            elif rest.mask1 and rest.mask2 and rest.mask3 and not rest.mask4:
                angle = pt.angle(traj, ' '.join([rest.mask1, rest.mask2, rest.mask3]))
                rest_dat[phase][-1].append(angle)
            elif rest.mask1 and rest.mask2 and rest.mask3 and rest.mask4:
                dihedral = pt.dihedral(traj, ' '.join([rest.mask1, rest.mask2, rest.mask3, rest.mask4]))
                rest_dat[phase][-1].append(dihedral)
            else:
                pass

    fecalc = paprika.analysis.fe_calc()
    fecalc.restraint_list = [rest1,rest2,rest3]
    fecalc.raw_values = rest_dat
    fecalc.compute_free_energy()

    # Test free energies and uncertainties
    test_vals = [
                fecalc.results['attach']['mbar']['blocking']['fe'],
                fecalc.results['attach']['mbar']['blocking']['sem'],
                fecalc.results['pull']['mbar']['blocking']['fe'],
                fecalc.results['pull']['mbar']['blocking']['sem']
                ]
    ref_vals = [13.267731176, 0.16892084090, -2.1791430735, 0.93638948302]

    assert np.allclose(test_vals, ref_vals)

    # Test convergence values attach

    test_vals = fecalc.results['attach']['mbar']['blocking']['convergence']

    ref_vals = np.array([0.0198918, 0.0451676, 0.0564517, 0.1079282, 0.1079282])

    assert np.allclose(test_vals, ref_vals)

    # Test convergence values pull

    test_vals = fecalc.results['pull']['mbar']['blocking']['convergence']

    ref_vals = np.array([0.2053769, 0.2053769, 0.1617423, 0.1747668, 0.5255023, 0.5255023, 0.1149945, 0.1707901, 0.2129136, 0.2129136, 0.1942189, 0.1768906, 0.1997338, 0.1997338, 0.2014766, 0.2014766, 0.1470727, 0.1442517, 0.1434395])

    assert np.allclose(test_vals, ref_vals)

#test_fe_calc()
