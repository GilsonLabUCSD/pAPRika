import parmed as pmd
import logging as log
import numpy as np

def align(structure, mask1, mask2, save=False, filename=None):
    """
    Align the mask1 -- mask2 vector to the z axis.
    """

    mask1_coordinates = structure[mask1].coordinates
    mask1_masses = [atom.mass for atom in structure[mask1].atoms]
    mask1_com = pmd.geometry.center_of_mass(np.asarray(mask1_coordinates),
                                            np.asarray(mask1_masses))

    mask2_coordinates = structure[mask2].coordinates
    mask2_masses = [atom.mass for atom in structure[mask2].atoms]
    mask2_com = pmd.geometry.center_of_mass(np.asarray(mask2_coordinates),
                                            np.asarray(mask2_masses))

    log.info('Moving {} ({} atoms) to the origin...'.format(mask1, len(mask1_coordinates)))
    log.info('Aligning {} ({} atoms) with the z axis...'.format(mask2, len(mask2_coordinates)))

    z = np.array([0.0, 0.0, 1.0])
    I = np.identity(3)
    # https://math.stackexchange.com/questions/293116/rotating-one-3d-vector-to-another
    # 1. Align mask1 to the origin.
    mask2_com = mask2_com + -1.0 * mask1_com

    # 2. Find axis and angle using cross and dot products.
    x     = np.cross(mask2_com, z) / np.linalg.norm(np.cross(mask2_com, z))
    theta = np.arccos(np.dot(mask2_com, z)/(np.linalg.norm(mask2_com) * np.linalg.norm(z)))
    # 3. Find the rotation matrix
    A = np.array([[0,        -1.0*x[2],     x[1] ],
                  [x[2],      0,       -1.0*x[0] ],
                  [-1.0*x[1], x[0],           0  ]])

    rotation_matrix = I + np.dot(np.sin(theta), A) + np.dot((1.0 - np.cos(theta)), np.dot(A, A))

    # This is certainly not the fastest approach, but it is explicit.
    aligned_coords = np.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        # Re-apply the translation of mask1
        aligned_coords[atom] = structure.coordinates[atom] + -1.0 * mask1_com
        aligned_coords[atom] = np.dot(rotation_matrix, aligned_coords[atom])
    structure.coordinates = aligned_coords

    if save:
        if not filename:
            log.warning('Unable to save aligned coordinates (no filename provided)...')
        else:
            log.info('Saved aligned coordinates to {}'.format(filename))
            # This seems to write out HETATM in place of ATOM
            # We should offer the option of writing a mol2 file, directly.
            structure.write_pdb(filename)

    return structure
    
def check_coordinates(structure, mask):
    """ 
    Return the coordinates of an atom selection.
    """
    mask_coordinates = structure[mask].coordinates
    mask_masses = [atom.mass for atom in structure[mask].atoms]
    mask_com = pmd.geometry.center_of_mass(np.asarray(mask_coordinates),
                                            np.asarray(mask_masses))
    return mask_com

def offset_structure(structure, offset):
    """
    Return a structure whose atoms have been translated.
    """
    offset_coords = np.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        offset_coords[atom] = structure.coordinates[atom] + offset
    structure.coordinates = offset_coords
    log.info('Added offset of {} to atomic coordinates...'.format(offset))
    return structure

def return_structure(filename):
    """
    Return structure object from file name.
    """
    # `parmed` can read both PDBs and
    # .inpcrd/.prmtop files with the same function call.
    try:
        structure = pmd.load_file(filename)
        log.info('Loaded {}...'.format(filename))
    except:
        log.error('Unable to load file: {}'.format(filename))
    return structure
