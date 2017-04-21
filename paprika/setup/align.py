import numpy as np
import parmed as pmd

def align(mask1, mask2, prmtop='test/cb6-but/vac.topo', inpcrd='test/cb6-but/vac.crds'):
    """
    Align the mask1 -- mask2 vector to the z axis.
    """

    structure = pmd.load_file(prmtop, inpcrd)
    mask1_coordinates = structure[mask1].coordinates
    mask1_masses = [atom.mass for atom in structure[mask1].atoms]
    mask1_com = pmd.geometry.center_of_mass(np.asarray(mask1_coordinates),
                                            np.asarray(mask1_masses))

    mask2_coordinates = structure[mask2].coordinates
    mask2_masses = [atom.mass for atom in structure[mask2].atoms]
    mask2_com = pmd.geometry.center_of_mass(np.asarray(mask2_coordinates),
                                            np.asarray(mask2_masses))

    print('Aligning {} ({} atoms) and {} ({} atoms) to the z axis.'.format(mask1,
                                                                           len(mask1_coordinates),
                                                                           mask2,
                                                                           len(mask2_coordinates)))
    print('COM of {} = {}'.format(mask1, mask1_com))
    print('COM of {} = {}'.format(mask2, mask2_com))

    z = np.array([0.0, 0.0, 1.0])
    I = np.identity(3)
    # https://math.stackexchange.com/questions/293116/rotating-one-3d-vector-to-another
    # 1. Align mask1 to the origin.
    mask1_com = mask1_com + -1.0 * mask1_com
    mask2_com = mask2_com + -1.0 * mask1_com
    print('{} should be at the origin: {}'.format(mask1, mask1_com))

    # 2. Find axis and angle using cross and dot products.
    x     = np.cross(mask2_com, z) / np.linalg.norm(np.cross(mask2_com, z))
    theta = np.arccos(np.dot(mask2_com, z)/(np.linalg.norm(mask2_com) * np.linalg.norm(z)))
    # 3. Find the rotation matrix
    A = np.array([[0,        -1.0*x[2],     x[1] ],
                  [x[2],      0,       -1.0*x[0] ],
                  [-1.0*x[1], x[0],           0  ]])

    rotation_matrix = I + np.dot(np.sin(theta), A) + np.dot((1.0 - np.cos(theta)), np.dot(A,A))

    return rotation_matrix
    
# coords = np.empty_like(structure.coordinates)
# for atom in range(len(structure.atoms)):
#     coords[atom] = np.dot(R, structure.coordinates[atom])
# # FYI, this will write out HETATM records.
# # This is not working yet.
# structure.coordinates = coords
# structure.write_pdb('tmp.pdb')
