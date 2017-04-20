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

    a = mask1_com
    b = mask2_com
    z = np.array([0.0, 0.0, 1.0])
    I = np.identity(3)
    # https://math.stackexchange.com/questions/293116/rotating-one-3d-vector-to-another
    # 1. Align mask1 to the origin.
    inversion = a * -1.0
    a = a + inversion
    b = b + inversion
    # 2. Find axis and angle using cross and dot products.
    x     = np.cross(b, z) / np.linalg.norm(np.cross(b, z))
    theta = np.arccos(np.dot(b, z)/(np.linalg.norm(b) * np.linalg.norm(z)))
    # 3. Find the rotation matrix
    A = np.array([[0,        -1.0*x[2],    x[1] ], 
                [ x[2],     0,       -1.0*x[0]], 
                [-1.0*x[1], x[0],          0] ])

    R = I + np.dot(np.sin(theta), A) + np.dot((1.0 - np.cos(theta)), np.dot(A,A))

    print(np.dot(R, structure.coordinates))