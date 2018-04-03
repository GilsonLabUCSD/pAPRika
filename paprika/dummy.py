import os as os
import re as re
import subprocess as sp

import logging as log
import numpy as np
import parmed as pmd
from parmed.structure import Structure as ParmedStructureClass
from paprika import utils


def add_dummy(structure,
              atom_name='DUM',
              residue_name='DUM',
              mass=208.00,
              atomic_number=82,
              x=0.000,
              y=0.000,
              z=0.000,
              ):
    """Add a dummy atom at the specified coordinates to the end of a structure.

    Parameters:
    ----------
    structure : str or pmd.Structure
        The structure to be modified
    atom_name : {str}, optional
        The name of the dummy atom (the default is 'DUM')
    residue_name : {str}, optional
        The residue name of the dummy atom (the default is 'DUM')
    mass : {float}, optional
        The mass of the dummy atom (the default is 208.00)
    atomic_number : {int}, optional
        The element of the dummy atom (the default is 82)
    x : {float}, optional
        The x coordinate of the dummy atom (the default is 0.000)
    y : {float}, optional
        The y coordinate of the dummy atom (the default is 0.000)
    z : {float}, optional
        The z coordinate of the dummy atom (the default is 0.000)
    mol2 : bool
        Whether to write a `mol2` file for the dummy atom
    frcmod : bool
        Whether to write a `frcmod` file for the dummy atom

    Returns:
    -------
    structure : pmd.Structure
        The modified structure
    """

    if isinstance(structure, str):
        structure = utils.return_parmed_structure(structure)
    elif isinstance(structure, ParmedStructureClass):
        pass
    else:
        raise Exception('add_dummy does not support the type associated with structure: ' + type(structure))

    # Create an atom object
    dum = pmd.topologyobjects.Atom()
    dum.name = atom_name
    dum.mass = mass
    dum.atomic_number = atomic_number
    # This may be a problem if these coordinates are outside the periodic box dimensions and ParmEd does not recalculate the box vectors before saving `inpcrd`...
    dum.xx = x
    dum.xy = y
    dum.xz = z

    # Assume that the last atom in the structure has the highest atom index, so the new atom will be at the end.
    dum.number = structure.atoms[-1].number + 1
    # Assume that the last residue in the structure has the highest residue number, so the new atom will be at the end.
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


def write_dummy_frcmod(atom_type='Du', mass='208.00', path='./', filename='dummy.frcmod', filepath=None):
    """Write a `frcmod` file for dummy atoms.

    Parameters:
    ----------
    atom_type : {str}, optional
        The atom type of the dummy atom (the default is 'Du')
    mass : {str}, optional
        The mass of the dummy atom (the default is '208.00')
    path : {str}, optional
        The directory of the output file, if `filepath` is not specified (the default is './')
    filename : {str}, optional
        The name of the output file, if `filepath` is not specified (the default is 'dummy.frcmod')
    filepath : {str}, optional
        The full path (directory and file) of the output (the default is None, which means `path` and `filename` will be used)

    """

    if filepath is None:
        filepath = os.path.join(path, filename)

    with open(filepath, 'w') as f:
        f.write("""\
Parameters for dummy atom with type {0}
MASS
{0}     {1}

BOND

ANGLE

DIHE

IMPROPER

NONBON
  {0}       0.000     0.0000000
""".format(atom_type, mass))


def write_dummy_mol2(atom_name='DUM',
                     atom_type='Du',
                     residue_name='DUM',
                     path='./',
                     filename='dummy.mol2',
                     filepath=None):
    """Write a `mol2` file for dummy atoms.

    Parameters:
    ----------
    atom_name : {str}, optional
        The atom name of the dummy atoms (the default is 'DUM')
    atom_type : {str}, optional
        The atom type of the dummy atoms (the default is 'Du')
    residue_name : {str}, optional
        The residue name of the dummy atoms (the default is 'DUM')
    path : {str}, optional
        The directory of the output file, if `filepath` is not specified (the default is './')
    filename : {str}, optional
        The name of the output file, if `filepath` is not specified (the default is 'dummy.mol2')
    filepath : {str}, optional
        The full path (directory and file) of the output (the default is None, which means `path` and `filename` will be used)

    """

    if filepath is None:
        filepath = os.path.join(path, filename)

    with open(filepath, 'w') as f:
        f.write("""\
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
""".format(residue_name[0:3], atom_name, atom_type[0:2]))
