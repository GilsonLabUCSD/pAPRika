import os as os

import parmed as pmd
from parmed.structure import Structure as ParmedStructureClass

from paprika import utils
from paprika.utils import check_unit
from openff.units import unit


def add_dummy(
    structure,
    atom_name="DUM",
    residue_name="DUM",
    mass=208.00 * unit.dalton,
    atomic_number=82,
    x=0.000 * unit.angstrom,
    y=0.000 * unit.angstrom,
    z=0.000 * unit.angstrom,
):
    """Add a dummy atom at the specified coordinates to the end of a structure.

    Parameters
    ----------
    structure : str or :class:`parmed.Structure`
        The structure to be modified
    atom_name : str, optional, default='DUM'
        The name of the dummy atom.
    residue_name : str, optional, default='DUM'
        The residue name of the dummy atom.
    mass : float, optional, default=208.0
        The mass of the dummy atom.
    atomic_number : int, optional, default=82
        The element of the dummy atom.
    x : float, optional, default=0.0
        The x coordinate of the dummy atom.
    y : float, optional, default=0.0
        The y coordinate of the dummy atom.
    z : float, optional, default=0.0
        The z coordinate of the dummy atom.

    Returns
    -------
    structure : :class:`parmed.Structure`
        The modified structure with dummy atoms added.

    """

    mass = check_unit(mass, base_unit=unit.dalton)
    x = check_unit(x, base_unit=unit.angstrom)
    y = check_unit(y, base_unit=unit.angstrom)
    z = check_unit(z, base_unit=unit.angstrom)

    if isinstance(structure, str):
        structure = utils.return_parmed_structure(structure)
    elif isinstance(structure, ParmedStructureClass):
        pass
    else:
        raise Exception(
            "add_dummy does not support the type associated with structure: "
            + type(structure)
        )

    # Create an atom object
    dum = pmd.topologyobjects.Atom()
    dum.name = atom_name
    dum.mass = mass.to(unit.dalton).magnitude
    dum.atomic_number = atomic_number
    # This may be a problem if these coordinates are outside the periodic box
    # dimensions and ParmEd does not recalculate the box vectors before saving
    # `inpcrd`...
    dum.xx = x.to(unit.angstrom).magnitude
    dum.xy = y.to(unit.angstrom).magnitude
    dum.xz = z.to(unit.angstrom).magnitude

    # Assume that the last atom in the structure has the highest atom index,
    # so the new atom will be at the end.
    dum.number = structure.atoms[-1].number + 1
    # Assume that the last residue in the structure has the highest residue
    # number, so the new atom will be at the end.
    residue_num = structure.residues[-1].number + 1

    structure.add_atom(dum, residue_name, residue_num)

    # Make sure that atom.number get set properly.  When reading in prmtop/inpcrd
    # parmed doesn't seem to set atom.number for some reason.
    for i, atom in enumerate(structure.atoms):
        atom.number = structure.atoms[i].idx + 1

    # tleap will probably want TER cards in any PDBs we make, so enforce
    # that for both the dummy residue and the residue before it
    structure.residues[-2].ter = True
    structure.residues[-1].ter = True

    return structure


def write_dummy_frcmod(
    atom_type="Du", mass="208.00", path="./", filename="dummy.frcmod", filepath=None
):
    """Write a ``frcmod`` file for dummy atoms.

    Parameters
    ----------
    atom_type : str, optional, default='Du'
        The atom type of the dummy atom.
    mass : str, optional, default=208.00
        The mass of the dummy atom.
    path : str, optional, default='./
        The directory of the output file, if `filepath` is not specified.
    filename : str, optional, default='dummy.frcmod`
        The name of the output file, if `filepath` is not specified.
    filepath : str, optional, default=None
        The full path (directory and file) of the output.

    """

    if filepath is None:
        filepath = os.path.join(path, filename)

    with open(filepath, "w") as f:
        f.write(
            """\
Parameters for dummy atom with type {0}
MASS
{0}     {1}

BOND

ANGLE

DIHE

IMPROPER

NONBON
  {0}       0.000     0.0000000
""".format(
                atom_type, mass
            )
        )


def write_dummy_mol2(
    atom_name="DUM",
    atom_type="Du",
    residue_name="DUM",
    path="./",
    filename="dummy.mol2",
    filepath=None,
):
    """Write a ``mol2`` file for dummy atoms.

    Parameters
    ----------
    atom_name : str, optional, default='DUM'
        The atom name of the dummy atoms.
    atom_type : str, optional, default='Du'
        The atom type of the dummy atoms.
    residue_name : str, optional, default='DUM'
        The residue name of the dummy atoms.
    path : str, optional, default='./'
        The directory of the output file, if `filepath` is not specified.
    filename : str, optional, default='dummy.mol2'
        The name of the output file, if `filepath` is not specified.
    filepath : str, optional, default=None
        The full path (directory and file) of the output.

    """

    if filepath is None:
        filepath = os.path.join(path, filename)

    with open(filepath, "w") as f:
        f.write(
            """\
@<TRIPOS>MOLECULE
{0}
    1     0     1     0     1
SMALL
USER_CHARGES

@<TRIPOS>ATOM
  1 {1:4s}    0.000000    0.000000    0.000000 {2}    1 {0}     0.0000 ****
@<TRIPOS>BOND
@<TRIPOS>SUBSTRUCTURE
      1  {0}              1 ****               0 ****  ****    0 ROOT
""".format(
                residue_name[0:3], atom_name, atom_type[0:2]
            )
        )


def extract_dummy_atoms(structure, resname=None, serial=True):
    """
    Extract information about dummy atoms from a parmed structure and
    returns the information as a dictionary.

    Parameters
    ----------
    structure : :class:`parmed.Structure`
        The parmed structure object we want to extract from
    resname : list
        List of residue name for the dummy atoms (default: ["DM1", "DM2", "DM3"])
    serial : bool
        get indices in serial (starts from 1) or index (starts from 0)

    Returns
    -------
    dummy_atoms : dict
        A dictionary containing positions (``pos``), index (``idx``) and mass (``mass``) of dummy atoms.
    """

    if resname is None:
        resname = ["DM1", "DM2", "DM3"]

    dummy_atoms = {name: {} for name in resname}

    for dummy_atom in resname:
        residue = f":{dummy_atom}"
        dummy_atoms[dummy_atom]["pos"] = structure[residue].coordinates[0]
        dummy_atoms[dummy_atom]["mass"] = [
            atom.mass for atom in structure[residue].atoms
        ][0]
        dummy_atoms[dummy_atom]["idx"] = utils.index_from_mask(
            structure, residue, amber_index=serial
        )[0]
        dummy_atoms[dummy_atom]["idx_type"] = "serial" if serial else "index"

    return dummy_atoms
