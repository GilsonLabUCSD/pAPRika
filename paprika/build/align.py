import logging

import numpy
import parmed
from openff.units import unit as openff_unit
from paprika.utils import check_unit

logger = logging.getLogger(__name__)


def zalign(structure, mask1, mask2, axis="z", save=False, filename=None):
    """Aligns the vector formed by atom mask1--mask2 to a reference axis.

    Parameters
    ----------
    structure : :class:`parmed.Structure`
        Molecular structure containing coordinates.
    mask1 : str
        Selection of first set of atoms.
    mask2 : str
        Selection of second set of atoms.
    axis : str or list or :class:`numpy.ndarray`, optional, default='z'
        Cartesian axis as a string: `x`, `y`, and `z`, or an arbitrary vector as a list. If an array is specified, the
        axis-vector will be normalized automatically.
    save : bool, optional, default=False
        Whether to save the coordinates (the default is False, which does nothing).
    filename : str, optional, default=None
        The filename for the saved coordinates (the default is None, which does nothing).

    Returns
    -------
    structure : :class:`parmed.Structure`
        A molecular structure with the coordinates aligned as specified.

    Examples
    --------
        >>> # Align a system based on the vector joining BUT@C--BUT@C3 to the z-axis
        >>> aligned_structure = zalign(structure, ":BUT@C", ":BUT@C3")
        >>>
        >>> # Do the same thing but align the system to the y-axis
        >>> aligned_structure = zalign(structure, ":BUT@C", ":BUT@C3", axis='y')
    """

    # Check axis
    axis = _return_array(axis)

    # Mask1
    mask1_coordinates = structure[mask1].coordinates
    mask1_masses = [atom.mass for atom in structure[mask1].atoms]
    mask1_com = parmed.geometry.center_of_mass(
        numpy.asarray(mask1_coordinates), numpy.asarray(mask1_masses)
    )

    # Mask2
    mask2_coordinates = structure[mask2].coordinates
    mask2_masses = [atom.mass for atom in structure[mask2].atoms]
    mask2_com = parmed.geometry.center_of_mass(
        numpy.asarray(mask2_coordinates), numpy.asarray(mask2_masses)
    )

    logger.info(
        "Moving {} ({} atoms) to the origin...".format(mask1, len(mask1_coordinates))
    )
    logger.info(
        "Aligning {} ({} atoms) with the z axis...".format(
            mask2, len(mask2_coordinates)
        )
    )

    # Define the vector from mask1 to mask2.
    mask2_com = mask2_com + -1.0 * mask1_com

    # Rotation matrix
    rotation_matrix = get_rotation_matrix(mask2_com, axis)

    # This is certainly not the fastest approach, but it is explicit.
    aligned_coords = numpy.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        aligned_coords[atom] = structure.coordinates[atom] + -1.0 * mask1_com
        aligned_coords[atom] = numpy.dot(rotation_matrix, aligned_coords[atom])
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


def align_principal_axes(structure, atom_mask=None, principal_axis=1, axis="z"):
    """Aligns a chosen principal axis of a system to a user specified axis. This function is based on the method given
    in the link: https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/

    Parameters
    ----------
    structure : :class:`parmed.Structure`
        Molecular structure containing coordinates.
    atom_mask : str, optional, default=None
        A mask that filter specific atoms for calculating the moment of inertia. This is useful if the molecule is
        not radially symmetric and you may want to select a subset of atoms that is symmetric for alignment.
    principal_axis : int, optional, default=1
        The particular principal axis to align to (The choices are `1`, `2` or `3` with `1` being the
        principal axis with the largest eigenvalue and `3` the lowest.).
    axis: str or list or :class:`numpy.ndarray`, optional, default='z'
        The axis vector to align the system to (by default the function aligns the principal axes with the
        largest eigenvalue to the z-axis). If an array is specified, the axis-vector will be normalized automatically.

    Returns
    -------
    structure : :class:`parmed.Structure`
        A molecular structure with it's principal axis aligned to a vector.

    Examples
    --------
        The commands below mimics the example given in the link above in VMD.

        >>> # Align largest principal axes to the z-axis
        >>> structure = align_principal_axes(structure, principal_axis=1, axis='z')
        >>>
        >>> # Align second largest principal axes to the y-axis
        >>> structure = align_principal_axes(structure, principal_axis=2, axis='y')

        It is also possible to use only a subset of atoms to determine the principal axes

        >>> structure = align_principal_axes(structure, atom_mask="@/C", principal_axis=2, axis='z')

    """
    # Check principal_axis
    if principal_axis is None:
        principal_axis = 1
    elif principal_axis not in [1, 2, 3]:
        raise ValueError(
            'Error: "principal_axis" can only be an integer value of 1, 2 or 3.'
        )

    # Check axis
    axis = _return_array(axis)

    # Principal axis vector
    p_axis = get_principal_axis_vector(structure, principal_axis, atom_mask)

    # Calculate Rotation matrix
    rotation_matrix = get_rotation_matrix(p_axis, axis)

    # Align the principal axis to specified axis
    aligned_coords = numpy.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        aligned_coords[atom] = structure.coordinates[atom]
        aligned_coords[atom] = numpy.dot(rotation_matrix, aligned_coords[atom])
    structure.coordinates = aligned_coords

    return structure


def rotate_around_axis(structure, axis, angle):
    """Rotates a system around an axis. The axis can be Cartesian axes ('x', 'y', and 'z') or an arbitrary axis
    (axis-angle).

    Parameters
    ----------
    structure : :class:`parmed.Structure`
        Molecular structure containing coordinates.
    axis: str or list or :class:`numpy.ndarray`
        The axis of rotation. If an array is specified, the axis-vector will be normalized automatically.
    angle: float or pint.unit.Quantity
        The angle of rotation in degrees.

    Returns
    -------
    structure : :class:`parmed.Structure`
        A molecular structure with it's principal axis aligned to a vector.

    Examples
    --------
        >>> # Rotate a system around the z axis by 90 degrees
        >>> structure = rotate_around_axis(structure, axis='z', angle=90.0)
        >>>
        >>> # Rotate a system around the axis [1,1,1] by 30 degrees
        >>> structure = rotate_around_axis(structure, axis=[1,1,1], angle=30.0)
    """

    # Check axis
    axis = _return_array(axis)

    # Convert angle to radians (temporary until Pint integration)
    angle = check_unit(angle, base_unit=openff_unit.radians).to(openff_unit.radians).magnitude

    if numpy.array_equal(axis, numpy.array([1.0, 0.0, 0.0])):
        rotation_matrix = numpy.array(
            [
                [1, 0, 0],
                [0, numpy.cos(angle), -numpy.sin(angle)],
                [0, numpy.sin(angle), numpy.cos(angle)],
            ]
        )
    elif numpy.array_equal(axis, numpy.array([0.0, 1.0, 0.0])):
        rotation_matrix = numpy.array(
            [
                [numpy.cos(angle), 0, numpy.sin(angle)],
                [0, 1, 0],
                [-numpy.sin(angle), 0, numpy.cos(angle)],
            ]
        )
    elif numpy.array_equal(axis, numpy.array([0.0, 0.0, 1.0])):
        rotation_matrix = numpy.array(
            [
                [numpy.cos(angle), -numpy.sin(angle), 0],
                [numpy.sin(angle), numpy.cos(angle), 0],
                [0, 0, 1],
            ]
        )
    else:
        # Normalize axis vector
        axis = axis / numpy.linalg.norm(axis)

        u_x = axis[0]
        u_y = axis[1]
        u_z = axis[2]

        rotation_matrix = numpy.array(
            [
                [
                    numpy.cos(angle) + u_x ** 2 * (1 - numpy.cos(angle)),
                    u_x * u_y * (1 - numpy.cos(angle)) - u_z * numpy.sin(angle),
                    u_x * u_z * (1 - numpy.cos(angle)) + u_y * numpy.sin(angle),
                ],
                [
                    u_y * u_x * (1 - numpy.cos(angle)) + u_z * numpy.sin(angle),
                    numpy.cos(angle) + u_y ** 2 * (1 - numpy.cos(angle)),
                    u_y * u_z * (1 - numpy.cos(angle)) - u_x * numpy.sin(angle),
                ],
                [
                    u_z * u_x * (1 - numpy.cos(angle)) - u_y * numpy.sin(angle),
                    u_z * u_y * (1 - numpy.cos(angle)) + u_x * numpy.sin(angle),
                    numpy.cos(angle) + u_z ** 2 * (1 - numpy.cos(angle)),
                ],
            ]
        )

    # Align the principal axis to specified axis
    aligned_coords = numpy.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        aligned_coords[atom] = structure.coordinates[atom]
        aligned_coords[atom] = numpy.dot(rotation_matrix, aligned_coords[atom])
    structure.coordinates = aligned_coords

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
    axis : str or list or :class:`numpy.ndarray`
        Cartesian axis as a string: `x`, `y`, and `z`, or an arbitrary vector as a list. If an array is specified,
        the axis vector need not be normalized.

    Returns
    -------
    theta : float
        The angle between the masks and the axis.
    """

    # Check axis
    axis = _return_array(axis)

    # Centroid of mask1
    mask1_coordinates = structure[mask1].coordinates
    mask1_masses = [atom.mass for atom in structure[mask1].atoms]
    mask1_com = parmed.geometry.center_of_mass(
        numpy.asarray(mask1_coordinates), numpy.asarray(mask1_masses)
    )

    # Centroid of mask2
    mask2_coordinates = structure[mask2].coordinates
    mask2_masses = [atom.mass for atom in structure[mask2].atoms]
    mask2_com = parmed.geometry.center_of_mass(
        numpy.asarray(mask2_coordinates), numpy.asarray(mask2_masses)
    )

    # Vector mask1-mask2
    vector = mask2_com + -1.0 * mask1_com

    # Angle between vectors
    theta = numpy.arccos(
        numpy.dot(vector, axis) / (numpy.linalg.norm(vector) * numpy.linalg.norm(axis))
    )

    return openff_unit.Quantity(theta, units=openff_unit.radians)


def get_rotation_matrix(vector, ref_vector):
    """Generate a rotation matrix needed for rotating ``vector`` to ``ref_vector``.

    Parameters
    ----------
    vector: :class:`numpy.ndarray`
        The vector of interest that will be rotated.
    ref_vector: :class:`numpy.ndarray`
        The reference vector, which vector1 will be rotated to.

    Returns
    -------
    rotation_matrix: :class:`numpy.ndarray`
        A 3x3 rotation matrix.

    """

    # If the structures are already aligned (cross product is zero), return 3x3 identity matrix
    if numpy.linalg.norm(numpy.cross(vector, ref_vector)) == 0:
        logger.info("The structure is already aligned and the denominator is invalid, returning identity matrix.")
        return numpy.identity(3)

    # Find axis between the mask vector and the axis using cross and dot products.
    x = numpy.cross(vector, ref_vector) / numpy.linalg.norm(numpy.cross(vector, ref_vector))

    theta = numpy.arccos(
        numpy.dot(vector, ref_vector)
        / (numpy.linalg.norm(vector) * numpy.linalg.norm(ref_vector))
    )

    # https://math.stackexchange.com/questions/293116/rotating-one-3d-vector-to-another
    A = numpy.array(
        [[0, -1.0 * x[2], x[1]], [x[2], 0, -1.0 * x[0]], [-1.0 * x[1], x[0], 0]]
    )

    # Rotation matrix
    rotation_matrix = (
        numpy.identity(3)
        + numpy.dot(numpy.sin(theta), A)
        + numpy.dot((1.0 - numpy.cos(theta)), numpy.dot(A, A))
    )

    return rotation_matrix


def get_principal_axis_vector(structure, principal_axis=1, atom_mask=None):
    """Return the principal axis vector given a structure.

    Parameters
    ----------
    structure: :class:`parmed.Structure`
        Molecular structure containing coordinates.
    principal_axis: int, optional, default=1
        The particular principal axis to align to (The choices are `1`, `2` or `3` with `1` being the
        principal axis with the largest eigenvalue and `3` the lowest).
    atom_mask: str, optional, default=None
        A mask that filter specific atoms for calculating the moment of inertia. This is useful if the molecule is
        not radially symmetric and you may want to select a subset of atoms that is symmetric for alignment.

    Returns
    -------
    p_axis: :class:`numpy.ndarray`
        Principal axis vector.
    """
    # Get coordinates and masses
    coordinates = structure.coordinates
    masses = numpy.asarray([atom.mass for atom in structure.atoms])
    if atom_mask:
        coordinates = structure[atom_mask].coordinates
        masses = numpy.asarray([atom.mass for atom in structure[atom_mask].atoms])

    # Calculate center of mass
    centroid = parmed.geometry.center_of_mass(coordinates, masses)

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

    inertia = numpy.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

    # Principal axis
    evals, evecs = numpy.linalg.eig(inertia)
    evecs = evecs[
        :, evals.argsort()[::-1]
    ]  # <-- Numpy may not sort the eigenvalues properly
    p_axis = evecs[:, principal_axis - 1]

    return p_axis


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
    mask_com : :class:`numpy.array`
        Coordinates of the selection center of mass.

    """

    mask_coordinates = structure[mask].coordinates
    mask_masses = [atom.mass for atom in structure[mask].atoms]
    mask_com = parmed.geometry.center_of_mass(
        numpy.asarray(mask_coordinates), numpy.asarray(mask_masses)
    )
    return mask_com


def offset_structure(structure, offset, dimension=None):
    """Return a structure whose coordinates have been shifted by ``offset``.

    Parameters
    ----------
    structure : :class:`parmed.Structure`
        Molecular structure containing coordinates.
    offset : float or :class:`numpy.ndarray` or pint.unit.Quantity
        The offset that will be added to *every* atom in the structure.
    dimension : str or list or :class:`numpy.ndarray`, optional, default=None
        By default the structure will be moved by ``offset`` in all direction if it is a single number, i.e. ``xyz +
        offset*[1,1,1]``. This variable applies a mask so that you can offset the structure in a specific direction.

    Returns
    -------
    structure : :class:`parmed.Structure`
        Coordinates of the structure offset by the given amount.

    Examples
    --------
        >>> # Offset a structure in the y-axis by 5 Angstrom
        >>> structure = offset_structure(structure, 5.0, dimension='y')
        >>>
        >>> # Offset a structure in the x- and z-axis by 3 Angstrom
        >>> structure = offset_structure(structure, 3.0, dimension='z')
        >>>
        >>> # Offset a structure using a numpy array
        >>> offset = numpy.array([0.0, 5.0, 2.0])
        >>> structure = offset_structure(structure, offset)
    """

    # Dimension mask
    if dimension is None:
        mask = numpy.array([1, 1, 1])
    else:
        mask = _return_array(dimension)

    offset = check_unit(offset, base_unit=openff_unit.angstrom)
    offset *= mask
    offset = offset.to(openff_unit.angstrom).magnitude

    # Offset coordinates
    offset_coords = numpy.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        offset_coords[atom] = structure.coordinates[atom] + offset
    structure.coordinates = offset_coords
    logger.info("Added offset of {} to atomic coordinates...".format(offset))

    return structure


def translate_to_origin(structure, weight="mass", atom_mask=None, dimension=None):
    """Translate a structure to the origin based on the centroid of the whole system or a subset of atom(s).

    Parameters
    ----------
    structure : str or :class:`parmed.Structure`
        Molecular structure containing coordinates.
    weight : str, optional, default="mass"
        Calculate the centroid based on either atomic masses (``mass`` default) or geometric center (``geo``).
    atom_mask : str, optional, default=None
        Selection of atom(s) if a particular subset is preferred to estimate the centroid.
    dimension : str or list or :class:`numpy.ndarray`, optional, default=None
        A mask that will filter the dimensions to which the translation will be applied (by default the system will
        be translated in all dimensions).

    Returns
    -------
    structure : :class:`parmed.Structure`
        A molecular structure with the coordinates translated to the origin.

    Examples
    --------
        >>> # Translate a structure to the origin
        >>> translated_structure = translate_to_origin(structure)
        >>>
        >>> # Translate the structure only in the z-axis
        >>> translated_structure = translate_to_origin(structure, dimension='z')
        >>>
        >>> # Translate the system based only on the carbon atoms
        >>> translated_structure = translate_to_origin(structure, atom_mask='@/C')
    """
    # Check if weight variable is properly chosen
    if weight not in ["mass", "geo"]:
        raise ValueError('"weight" must either be "mass" or "geo".')

    # Dimension mask
    if dimension is None:
        mask = numpy.array([1, 1, 1])
    else:
        mask = _return_array(dimension)

    # Atomic coordinates and masses
    if atom_mask is None:
        coordinates = structure.coordinates
        masses = numpy.asarray([atom.mass for atom in structure.atoms])
    else:
        coordinates = structure[atom_mask].coordinates
        masses = numpy.asarray([atom.mass for atom in structure[atom_mask].atoms])

    # Equal weights if geometric center is preferred
    if weight == "geo":
        masses = numpy.ones(len(coordinates))

    # Centroid coordinates
    centroid = parmed.geometry.center_of_mass(coordinates, masses)

    if mask is not None:
        centroid *= mask

    # Translate coordinates
    aligned_coords = numpy.empty_like(structure.coordinates)
    for atom in range(len(structure.atoms)):
        aligned_coords[atom] = structure.coordinates[atom] - centroid
    structure.coordinates = aligned_coords

    return structure


def _return_array(var):
    """Private module to return an array given a string, list or numpy array."""

    array = None

    if isinstance(var, str):
        if "x" in var.lower():
            array = numpy.array([1.0, 0.0, 0.0])
        elif "y" in var.lower():
            array = numpy.array([0.0, 1.0, 0.0])
        elif "z" in var.lower():
            array = numpy.array([0.0, 0.0, 1.0])
        else:
            raise KeyError(f"Cannot understand the variable: {var}.")

    elif isinstance(var, list):
        if len(var) != 3:
            raise ValueError('"var" must be a list with 3 elements.')
        array = numpy.array(var)

    elif isinstance(var, numpy.ndarray):
        if len(var) != 3:
            raise ValueError('"var" must be a numpy array with 3 elements.')

    else:
        raise TypeError(f"Variable type {type(var)} is not supported.")

    return array
