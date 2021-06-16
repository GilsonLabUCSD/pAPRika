import sys

import click
import parmed as pmd

from .utils import OrderedGroup


# Commands for `Align` modules
@click.group(cls=OrderedGroup)
def align():
    """Command-line interface for the `Align` modules."""


@align.command()
@click.option(
    "-c",
    "--coordinates",
    "coordinates_file",
    required=True,
    type=str,
    help="The coordinates file (pdb, rst7, mol2, gro)",
)
@click.option(
    "-t",
    "--topology",
    "topology_file",
    default=None,
    type=str,
    help="The topology file (pdb, psf, prmtop, top).",
)
@click.option(
    "-o",
    "--output",
    "output_file",
    default="output",
    type=str,
    help="The name for the output file (default: output.pdb).",
)
@click.option(
    "-m1",
    "--mask1",
    "mask1",
    required=True,
    help="Selection of first set of atoms.",
)
@click.option(
    "-m2",
    "--mask2",
    "mask2",
    required=True,
    help="Selection of second set of atoms.",
)
@click.option(
    "--axis",
    "axis",
    default="z",
    type=str,
    help="The reference Cartesian axis to align the structure (x, y, or z) (default: z).",
)
@click.option(
    "--vector",
    "vector",
    nargs=3,
    type=float,
    default=None,
    help="An arbitrary reference axis in vector form (e.g. 1 1 1).",
)
def zalign(
    coordinates_file,
    topology_file,
    output_file,
    mask1,
    mask2,
    axis,
    vector,
):
    """Aligns a structure to an axis based on vec(mask1--mask2)."""
    from paprika.build.align import zalign

    # Load files as a ParmEd structure
    files = [coordinates_file]
    if topology_file is not None:
        files = [topology_file] + files

    structure = pmd.load_file(*files, structure=True)

    # Translate structure
    if vector is not None:
        print(f"zalign: Aligning vec({mask1}-{mask2}) to the {axis}-axis.")
        structure = zalign(structure, mask1, mask2, axis=axis)
    else:
        print(f"zalign: Aligning vec({mask1}-{mask2}) to {vector}.")
        structure = zalign(structure, mask1, mask2, axis=list(vector))

    # Save to files
    structure.save(output_file, overwrite=True)


@align.command()
@click.option(
    "-c",
    "--coordinates",
    "coordinates_file",
    required=True,
    type=str,
    help="The coordinates file (pdb, rst7, mol2, gro)",
)
@click.option(
    "-t",
    "--topology",
    "topology_file",
    default=None,
    type=str,
    help="The topology file (pdb, psf, prmtop, top).",
)
@click.option(
    "-o",
    "--output",
    "output_file",
    default="output",
    type=str,
    help="The name for the output file (default: output.pdb).",
)
@click.option(
    "--p_axis",
    "p_axis",
    default=1,
    type=int,
    help="The choice of principal axis for the alignment (1, 2, or 3) (default: 1, i.e. largest axis).",
)
@click.option(
    "--axis",
    "axis",
    default="z",
    type=str,
    help="The reference Cartesian axis to align the structure (x, y, or z) (default: z).",
)
@click.option(
    "--vector",
    "vector",
    nargs=3,
    type=float,
    default=None,
    help="An arbitrary reference axis in vector form (e.g. 1 1 1).",
)
@click.option(
    "--mask",
    "atom_mask",
    type=str,
    default=None,
    help="Selection of atoms for calculating the moment of inertia (default: all atoms).",
)
def princ(
    coordinates_file,
    topology_file,
    output_file,
    p_axis,
    axis,
    vector,
    atom_mask,
):
    """Aligns a principal-axis of the structure to an reference axis."""
    from paprika.build.align import align_principal_axes

    # Load files as a ParmEd structure
    files = [coordinates_file]
    if topology_file is not None:
        files = [topology_file] + files

    structure = pmd.load_file(*files, structure=True)

    # Align structure
    if len(vector) != 3:
        print(f"princ: Aligning principal axis [{p_axis}] to the {axis}-axis.")
        structure = align_principal_axes(
            structure, atom_mask=atom_mask, principal_axis=p_axis, axis=axis
        )
    else:
        print(f"princ: Aligning principal axis [{p_axis}] to {vector}.")
        structure = align_principal_axes(
            structure, atom_mask=atom_mask, principal_axis=p_axis, axis=list(vector)
        )

    # Save to files
    structure.save(output_file, overwrite=True)


@align.command()
@click.option(
    "-c",
    "--coordinates",
    "coordinates_file",
    required=True,
    type=str,
    help="The coordinates file (pdb, rst7, mol2, gro)",
)
@click.option(
    "-t",
    "--topology",
    "topology_file",
    default=None,
    type=str,
    help="The topology file (pdb, psf, prmtop, top).",
)
@click.option(
    "-o",
    "--output",
    "output_file",
    required=True,
    type=str,
    help="The name for the output file (default: output.pdb).",
)
@click.option(
    "--mask",
    "atom_mask",
    default=None,
    help="The selection of atoms as reference when translating the structure to the origin (default: all atoms).",
)
@click.option(
    "--origin", "origin", is_flag=True, help="Translate the structure to the origin."
)
@click.option(
    "--offset",
    "offset",
    default=None,
    type=float,
    help="Offset the structure by this amount (direction needs to be specified).",
)
@click.option(
    "--direction",
    "direction",
    default=None,
    type=str,
    help="The direction to offset the structure +-(x, y, or z).",
)
@click.option(
    "--vector",
    "vector",
    nargs=3,
    type=float,
    default=None,
    help="Translate the structure with a vector.",
)
def translate(
    coordinates_file,
    topology_file,
    output_file,
    atom_mask,
    origin,
    offset,
    direction,
    vector,
):
    """Translate a structure to the origin or by an arbitrary amount."""
    import numpy as np

    from paprika.build.align import offset_structure, translate_to_origin

    # Load files as a ParmEd structure
    files = [coordinates_file]
    if topology_file is not None:
        files = [topology_file] + files

    structure = pmd.load_file(*files, structure=True)

    # Translate structure
    if offset is not None and direction is None:
        print("error: `axis` needs to be defined if `offset` is specified.")
        sys.exit()
    elif offset is not None and direction is not None:
        print(
            f"translate: Translating structure by {offset} in the `{direction}` direction."
        )
        structure = offset_structure(structure, offset, dimension=direction)
    if origin:
        print("translate: Translating structure to the origin.")
        structure = translate_to_origin(structure, atom_mask=atom_mask, dimension=direction)
    elif not origin and len(vector) == 3:
        print(f"translate: Translating structure by {vector}.")
        structure = offset_structure(structure, np.array(vector), dimension=None)

    # Save to files
    structure.save(output_file, overwrite=True)


@align.command()
@click.option(
    "-c",
    "--coordinates",
    "coordinates_file",
    required=True,
    type=str,
    help="The coordinates file (pdb, rst7, mol2, gro)",
)
@click.option(
    "-t",
    "--topology",
    "topology_file",
    default=None,
    type=str,
    help="The topology file (pdb, psf, prmtop, top).",
)
@click.option(
    "-o",
    "--output",
    "output_file",
    default="output",
    type=str,
    help="The name for the output file (default: output.pdb).",
)
@click.option(
    "-a",
    "--angle",
    "angle",
    required=True,
    type=float,
    help="The angle of rotation in degrees.",
)
@click.option(
    "-x",
    "--axis",
    "axis",
    type=str,
    default=None,
    help="The axis of rotation (x, y, or z).",
)
@click.option(
    "-v",
    "--vector",
    "vector",
    nargs=3,
    type=float,
    default=None,
    help="The axis of rotation in vector form (e.g. 1 1 1). "
    "This option supersedes the choice of axis if both `axis` and `vector` are specified.",
)
def rotate(
    coordinates_file, topology_file, output_file, angle, axis, vector
):
    """Rotate a structure around an axis."""
    from paprika.build.align import rotate_around_axis

    # Load files as a ParmEd structure
    files = [coordinates_file]
    if topology_file is not None:
        files = [topology_file] + files

    structure = pmd.load_file(*files, structure=True)

    # Rotate structure
    if axis is not None and len(vector) != 3:
        print(f"rotate: Rotating structure around the {axis}-axis by {angle} degrees.")
        structure = rotate_around_axis(structure, axis=axis, angle=angle)
    elif axis is None and len(vector) == 3:
        print(f"rotate: Rotating structure around {vector} by {angle} degrees.")
        structure = rotate_around_axis(structure, axis=list(vector), angle=angle)
    else:
        print(f"rotate: Rotating structure around {vector} by {angle} degrees.")
        structure = rotate_around_axis(structure, axis=list(vector), angle=angle)

    # Save to file
    structure.save(output_file, overwrite=True)


# Commands for `Dummy Atom` modules
@click.group("dummy", cls=OrderedGroup)
def dummy_atoms():
    """Command-line interface for the `Dummy Atom` modules."""


@dummy_atoms.command()
@click.option(
    "-c",
    "--coordinates",
    "coordinates_file",
    required=True,
    type=str,
    help="The coordinates file (pdb, rst7, mol2, gro)",
)
@click.option(
    "-t",
    "--topology",
    "topology_file",
    default=None,
    type=str,
    help="The topology file (pdb, psf, prmtop, top).",
)
@click.option(
    "-o",
    "--output",
    "output_name",
    required=True,
    type=str,
    help="The prefix name for the output.",
)
@click.option(
    "-fo",
    "--file_format",
    "file_format",
    nargs=0,
    default=["pdb"],
    help="The file(s) to export (pdb, mol2, gro, rst7, top, psf, prmtop).",
)
@click.argument("file_format", nargs=-1)
@click.option(
    "-x",
    "--x_pos",
    "x_pos",
    default=0.0,
    type=float,
    help="The x position of the Dummy atom (default: 0.0).",
)
@click.option(
    "-y",
    "--y_pos",
    "y_pos",
    default=0.0,
    type=float,
    help="The y position of the Dummy atom (default: 0.0).",
)
@click.option(
    "-z",
    "--z_pos",
    "z_pos",
    default=0.0,
    type=float,
    help="The z position of the Dummy atom (default: 0.0).",
)
def add(coordinates_file, topology_file, output_name, file_format, x_pos, y_pos, z_pos):
    """Add dummy atoms to a structure file."""
    from paprika.build.dummy import add_dummy

    # Load files to ParmEd structure
    files = [coordinates_file]
    if topology_file is not None:
        files = [topology_file] + files

    structure = pmd.load_file(*files, structure=True)

    # Add dummy atoms
    structure = add_dummy(structure, residue_name="DM1", x=x_pos, y=y_pos, z=z_pos)

    # Save to files
    for ext in file_format:
        structure.save(f"{output_name}.{ext.lower()}", overwrite=True)


@dummy_atoms.command()
@click.option(
    "-o",
    "--output",
    "output",
    required=True,
    type=str,
    help="The name for the frcmod file.",
)
@click.option(
    "-m",
    "--mass",
    "mass",
    default=208.0,
    type=float,
    help="The mass of the Dummy atom (default: 208.0 daltons).",
)
@click.option(
    "-t",
    "--atom_type",
    "atom_type",
    default="Du",
    type=str,
    help="The atom type name of the Dummy atom (default: 'Du').",
)
def frcmod(atom_type, output, mass):
    """Write a frcmod file for Dummy atom."""
    from paprika.build.dummy import write_dummy_frcmod

    write_dummy_frcmod(atom_type=atom_type, mass=mass, filename=output)


@dummy_atoms.command()
@click.option(
    "-o",
    "--output",
    "output_name",
    required=True,
    type=str,
    help="The name for the MOL2 file.",
)
@click.option(
    "-t",
    "--atom_type",
    "atom_type",
    default="Du",
    type=str,
    help="The atom type name of the Dummy atom (default: 'Du').",
)
@click.option(
    "-n",
    "--atom_name",
    "atom_name",
    default="DUM",
    type=str,
    help="The atom name for the Dummy atom (default: 'DUM').",
)
@click.option(
    "-r",
    "--residue_name",
    "residue_name",
    default="DUM",
    type=str,
    help="The residue name for the Dummy atom (default: 'DUM').",
)
def mol2(atom_type, atom_name, residue_name, output_name):
    """Write a MOL2 file for Dummy atom."""
    from paprika.build.dummy import write_dummy_mol2

    write_dummy_mol2(
        atom_name=atom_name,
        atom_type=atom_type,
        residue_name=residue_name,
        path="./",
        filename=output_name,
        filepath=None,
    )
