import logging

import numpy as np
import parmed as pmd

logger = logging.getLogger(__name__)


def zalign(structure, mask1, mask2, save=False, filename=None):
    """Aligns the vector formed by atom mask1--mask2 to the z-axis.

    Parameters
    ----------
    structure : :class:`parmed.Structure`
        Molecular structure containing coordinates.
    mask1 : str
        Selection of first set of atoms.
    mask2 : str
        Selection of second set of atoms.
    save : bool, optional
        Whether to save the coordinates (the default is False, which does nothing).
    filename : str, optional
        The filename for the saved coordinates (the default is None, which does nothing).

    Returns
    -------
    structure : :class:`parmed.Structure`
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
    """Get the angle (theta) between the vector formed by atom mask1--mask2 and a Cartesian axis.

    Parameters
    ----------
    structure : :class:`parmed.Structure`
        Molecular structure containing coordinates.
    mask1 : str
        Selection of first set of atoms.
    mask2 : str
        Selection of second set of atoms.
    axis : str
        Axis: x, y, or z.

    Returns
    -------
    theta : float
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
    structure : :class:`parmed.Structure`
        Molecular structure containing coordinates.
    mask : str
        Amber-style atom selection.

    Returns
    -------
    mask_com : :class:`np.array`
        Coordinates of the selection center of mass.

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
    structure : :class:`parmed.Structure`
        Molecular structure containing coordinates.
    offset : float
        The offset that will be added to *every* atom in the structure.

    Returns
    -------
    structure : :class:`parmed.Structure`
        Coordinates of the structure offset by the given amount.

    """

    offset_coords = np.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        offset_coords[atom] = structure.coordinates[atom] + offset
    structure.coordinates = offset_coords
    logger.info("Added offset of {} to atomic coordinates...".format(offset))
    return structure


def translate_to_origin(structure, weight="mass", atom_mask=None, dimension_mask=None):
    """Translate a structure to the origin based on the centroid of the whole system or a subset of atom(s).

    Parameters
    ----------
    structure : str or :class:`parmed.Structure`
        Molecular structure containing coordinates.
    weight : str, optional, default="mass"
        Calculate the centroid based on either atomic masses (``mass`` default) or geometric center (``geo``).
    atom_mask : str, optional, default=None
        Selection of atom(s) if a particular subset is preferred to estimate the centroid.
    dimension_mask : list, optional, default=None
        A mask that will filter the dimensions to which the translation will be applied (by default the system will
        be translated in all dimensions). For example, ``dimension_mask=[0,0,1]`` will translate the system to the
        origin only in the z-axis.

    Returns
    -------
    structure : :class:`parmed.Structure`
        A molecular structure with the coordinates translated to the origin.

    """
    # Check if weight variable is properly chosen
    if weight not in ["mass", "geo"]:
        raise Exception('Error: "weight" must either be "mass" or "geo".')

    # Dimension mask
    if dimension_mask is None:
        dimension_mask = [1, 1, 1]
    elif len(dimension_mask) != 3:
        raise Exception(
            'Error: "dimension_mask" must be a list with 3 elements, e.g. [1,1,1].'
        )
    dimension_mask = np.array(dimension_mask)

    # Atomic coordinates and masses
    if atom_mask is None:
        coordinates = structure.coordinates
        masses = np.asarray([atom.mass for atom in structure.atoms])
    else:
        coordinates = structure[atom_mask].coordinates
        masses = np.asarray([atom.mass for atom in structure[atom_mask].atoms])

    # Equal weights if geometric center is preferred
    if weight == "geo":
        masses[:] = 1.0

    # Centroid coordinates
    centroid = pmd.geometry.center_of_mass(coordinates, masses)

    # Translate coordinates
    aligned_coords = np.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        aligned_coords[atom] = structure.coordinates[atom] - centroid * dimension_mask
    structure.coordinates = aligned_coords

    return structure


def align_principal_axes(structure, atom_mask=None, principal_axis=1, v_axis=None):
    """Aligns the a chosen principal axis of a system to a specified axis. This function is based on the method given
    in the link: https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/

    Parameters
    ----------
    structure : :class:`parmed.Structure`
        Molecular structure containing coordinates.
    atom_mask : str, optional, default=None
        A mask that filter specific atoms for calculating the moment of inertia.
    principal_axis : int, optional, default=1
        The particular principal axis to align to (The choices are `1`, `2` or `3` with `1` being the
        principal axis with the largest eigenvalue and `3` the lowest.).
    v_axis: list, optional, default=None
        The axis vector to align the system to (by default the function aligns the principal axes with the
        largest eigenvalue to the z-axis).

    Returns
    -------
    structure : :class:`parmed.Structure`
        A molecular structure with it's principal axis aligned to a vector.

    Examples
    --------
    The commands below mimics the example given in the link above in VMD.

        >>> structure = align_principal_axes(structure, princ_axis=0, v_axis=[0,0,1])
        >>> structure = align_principal_axes(structure, princ_axis=1, v_axis=[0,1,0])

    """
    # Check principal_axis
    if principal_axis is None:
        principal_axis = 1
    elif principal_axis not in [0, 1, 2]:
        raise Exception(
            'Error: "principal_axis" can only be an integer value of 1, 2 or 3.'
        )

    # Axis vector for alignment
    if v_axis is None:
        v_axis = [0, 0, 1]
    elif len(v_axis) != 3:
        raise Exception('Error: "v_axis" must be a list with 3 elements, e.g. [0,0,1]')
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
    p_axis = evecs[:, principal_axis]

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
