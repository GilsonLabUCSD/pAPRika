import logging

import numpy as np
import parmed as pmd

logger = logging.getLogger(__name__)


def zalign(structure, mask1, mask2, save=False, filename=None):
    """Align the mask1 -- mask2 vector to the z axis.

    Parameters
    ----------
    structure : parmed.Structure
        Molecular structure containing coordinates
    mask1 : str
        Selection of first set of atoms
    mask2 : str
        Selection of second set of atoms
    save : bool, optional
        Whether to save the coordinates (the default is False, which does nothing)
    filename : str, optional
        The filename for the saved coordinates (the default is None, which does nothing)

    Returns
    -------
    parmed.Structure
        A molecular structure with the coordinates aligned as specified.
    """

    mask1_coordinates = structure[mask1].coordinates
    mask1_masses = [atom.mass for atom in structure[mask1].atoms]
    mask1_com = pmd.geometry.center_of_mass(
        np.asarray(mask1_coordinates), np.asarray(mask1_masses)
    )

    mask2_coordinates = structure[mask2].coordinates
    mask2_masses = [atom.mass for atom in structure[mask2].atoms]
    mask2_com = pmd.geometry.center_of_mass(
        np.asarray(mask2_coordinates), np.asarray(mask2_masses)
    )

    logger.info(
        "Moving {} ({} atoms) to the origin...".format(mask1, len(mask1_coordinates))
    )
    logger.info(
        "Aligning {} ({} atoms) with the z axis...".format(
            mask2, len(mask2_coordinates)
        )
    )

    axis = np.array([0.0, 0.0, 1.0])

    identity = np.identity(3)
    # https://math.stackexchange.com/questions/293116/rotating-one-3d-vector-to-another
    # 1. Define the vector from mask1 to mask2.
    mask2_com = mask2_com + -1.0 * mask1_com

    # 2. Find axis and angle between the mask vector and the axis using cross
    # and dot products.
    try:
        x = np.cross(mask2_com, axis) / np.linalg.norm(np.cross(mask2_com, axis))
    except RuntimeWarning:
        # The structure is already aligned and the denominator is invalid
        pass

    theta = np.arccos(
        np.dot(mask2_com, axis) / (np.linalg.norm(mask2_com) * np.linalg.norm(axis))
    )
    # 3. Find the rotation matrix
    A = np.array(
        [[0, -1.0 * x[2], x[1]], [x[2], 0, -1.0 * x[0]], [-1.0 * x[1], x[0], 0]]
    )

    rotation_matrix = (
        identity
        + np.dot(np.sin(theta), A)
        + np.dot((1.0 - np.cos(theta)), np.dot(A, A))
    )

    # This is certainly not the fastest approach, but it is explicit.
    aligned_coords = np.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        aligned_coords[atom] = structure.coordinates[atom] + -1.0 * mask1_com
        aligned_coords[atom] = np.dot(rotation_matrix, aligned_coords[atom])
    structure.coordinates = aligned_coords

    if save:
        if not filename:
            logger.warning(
                "Unable to save aligned coordinates (no filename provided)..."
            )
        else:
            logger.info("Saved aligned coordinates to {}".format(filename))
            # This seems to write out HETATM in place of ATOM
            # We should offer the option of writing a mol2 file, directly.
            structure.write_pdb(filename)

    return structure


def get_theta(structure, mask1, mask2, axis):
    """Get the angle (theta) between the vector formed by two masks and an axis.
    
    Parameters
    ----------
    structure : parmed.Structure
        Molecular structure containing coordinates
    mask1 : str
        Selection of first set of atoms
    mask2 : str
        Selection of second set of atoms
    axis : str
        Axis: x, y, or z
    
    Returns
    -------
    float
        The angle between the masks and the axis.
    """

    if "x" in axis.lower():
        axis = np.array([1.0, 0.0, 0.0])
    elif "y" in axis.lower():
        axis = np.array([0.0, 1.0, 0.0])
    elif "z" in axis.lower():
        axis = np.array([0.0, 0.0, 1.0])

    mask1_coordinates = structure[mask1].coordinates
    mask1_masses = [atom.mass for atom in structure[mask1].atoms]
    mask1_com = pmd.geometry.center_of_mass(
        np.asarray(mask1_coordinates), np.asarray(mask1_masses)
    )

    mask2_coordinates = structure[mask2].coordinates
    mask2_masses = [atom.mass for atom in structure[mask2].atoms]
    mask2_com = pmd.geometry.center_of_mass(
        np.asarray(mask2_coordinates), np.asarray(mask2_masses)
    )

    vector = mask2_com + -1.0 * mask1_com
    theta = np.arccos(
        np.dot(vector, axis) / (np.linalg.norm(vector) * np.linalg.norm(axis))
    )

    return theta


def check_coordinates(structure, mask):
    """Return the coordinates of an atom selection.
    
    Parameters
    ----------
    structure : parmed.Structure
        Molecular structure containing coordinates
    mask : str
        Atom selection
    
    Returns
    -------
    np.array
        Coordinates of the selection center of mass
    """

    mask_coordinates = structure[mask].coordinates
    mask_masses = [atom.mass for atom in structure[mask].atoms]
    mask_com = pmd.geometry.center_of_mass(
        np.asarray(mask_coordinates), np.asarray(mask_masses)
    )
    return mask_com


def offset_structure(structure, offset):
    """Return a structure whose coordinates have been offset.
    
    Parameters
    ----------
    structure : parmed.Structure
        Molecular structure containing coordinates
    offset : float
        The offset that will be added to *every* atom in the structure
    
    Returns
    -------
    :py:class:`parmed.Structure`
        Coordinates of the structure offset by the given amount.
    """

    offset_coords = np.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        offset_coords[atom] = structure.coordinates[atom] + offset
    structure.coordinates = offset_coords
    logger.info("Added offset of {} to atomic coordinates...".format(offset))
    return structure


def translate_to_origin(structure, weight="mass", atom_mask=None, dim_mask=None):
    """Translate a structure to the origin based on the center of 
       mass of the whole system or of a particular choice of atom(s).

    Parameters
    ----------
    structure : str or parmed.Structure
        Molecular structure containing coordinates.
    weight : str
        Calculate the center based on atomic masses ('mass' default) or geometric center ('geo')?
    atom_mask : str
        Selection of atom(s) if a particular subset is preferred to estimate the center.
    dim_mask : list
        A mask that will filter the dimensions to which the translation will be applied (default is [1,1,1]).

    Returns
    -------
    structure : parmed.Structure
        A molecular structure with the coordinates translated to the origin.

    """
    # Check if weight variable is properly chosen
    if weight not in ["mass", "geo"]:
        raise SystemExit('Error: "weight" must either be "mass" or "geo"')

    # Dimension mask
    if dim_mask is None:
        dim_mask = [1, 1, 1]
    elif len(dim_mask) != 3:
        raise SystemExit(
            'Error: "dim_mask" must be a list with 3 elements, e.g. [1,1,1]'
        )
    dim_mask = np.array(dim_mask)

    # Atom coordinates and masses
    if atom_mask is None:
        coordinates = structure.coordinates
        masses = np.asarray([atom.mass for atom in structure.atoms])
    else:
        coordinates = structure[atom_mask].coordinates
        masses = np.asarray([atom.mass for atom in structure[atom_mask].atoms])

    # Convert weight if geometric center
    if weight == "geo":
        masses[:] = 1.0

    # Center of mass/geometry
    centroid = pmd.geometry.center_of_mass(coordinates, masses)

    # Translate coordinates
    aligned_coords = np.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        aligned_coords[atom] = structure.coordinates[atom] - centroid * dim_mask
    structure.coordinates = aligned_coords

    return structure


def align_principal_axes(structure, atom_mask=None, princ_axis=None, v_axis=None):
    """Aligns the a chosen principal axis of a system to a specified axis. 
       This function is based on the method given in the link:

        * https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/

    Parameters
    ----------
    structure : parmed.Structure
        Molecular structure containing coordinates.
    atom_mask : str
        A mask that filter specific atoms for calculating the moment of inertia.
    princ_axis : int
        Which principal axis to align? (choice: 1, 2 or 3, default: 1).
    v_axis: list
        The axis vector to align the system to (default: [0, 0, 1]).

    Returns
    -------
    structure : parmed.Structure
        A molecular structure with a principal axis aligned to a vector.

    Examples
    --------
    The commands below mimics the example given in the link above for VMD.

        >>> structure = align_principal_axes(structure, princ_axis=0, v_axis=[0,0,1])
        >>> structure = align_principal_axes(structure, princ_axis=1, v_axis=[0,1,0])

    """
    # Check princ_axis
    if princ_axis is None:
        princ_axis = 1
    elif princ_axis not in [0, 1, 2]:
        raise SystemExit(
            'Error: "princ_axis" can only be an integer value of 0, 1 or 2.'
        )

    # Axis vector for alignment
    if v_axis is None:
        v_axis = [0, 0, 1]
    elif len(v_axis) != 3:
        raise SystemExit('Error: "v_axis" must be a list with 3 elements, e.g. [0,0,1]')
    v_axis = np.array(v_axis)

    # Get coordinates and masses
    coordinates = structure.coordinates
    masses = np.asarray([atom.mass for atom in structure.atoms])
    if atom_mask:
        coordinates = structure[atom_mask].coordinates
        masses = np.asarray([atom.mass for atom in structure[atom_mask].atoms])

    # Calculate center of mass
    centroid = pmd.geometry.center_of_mass(coordinates, masses)

    # Construct Inertia tensor
    Ixx = Ixy = Ixz = Iyy = Iyz = Izz = 0
    for xyz, mass in zip(coordinates, masses):
        xyz -= centroid
        Ixx += mass * (xyz[1] * xyz[1] + xyz[2] * xyz[2])
        Ixy -= mass * (xyz[0] * xyz[1])
        Ixz -= mass * (xyz[0] * xyz[2])
        Iyy += mass * (xyz[0] * xyz[0] + xyz[2] * xyz[2])
        Iyz -= mass * (xyz[1] * xyz[2])
        Izz += mass * (xyz[0] * xyz[0] + xyz[1] * xyz[1])

    inertia = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

    # Principal axis
    evals, evecs = np.linalg.eig(inertia)
    evecs = evecs[:, evals.argsort()]  # <-- Numpy does not sort the vectors properly
    p_axis = evecs[:, princ_axis]

    # Calculate Rotation matrix
    x = np.cross(p_axis, v_axis) / np.linalg.norm(np.cross(p_axis, v_axis))
    theta = np.arccos(
        np.dot(p_axis, v_axis) / (np.linalg.norm(p_axis) * np.linalg.norm(v_axis))
    )
    A = np.array(
        [[0, -1.0 * x[2], x[1]], [x[2], 0, -1.0 * x[0]], [-1.0 * x[1], x[0], 0]]
    )
    rotation_matrix = (
        np.identity(3)
        + np.dot(np.sin(theta), A)
        + np.dot((1.0 - np.cos(theta)), np.dot(A, A))
    )

    # Align the principal axis to specified axis
    aligned_coords = np.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        aligned_coords[atom] = structure.coordinates[atom]
        aligned_coords[atom] = np.dot(rotation_matrix, aligned_coords[atom])
    structure.coordinates = aligned_coords

    return structure
