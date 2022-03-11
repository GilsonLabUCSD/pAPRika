from enum import Enum

N_A = 6.0221409 * 10 ** 23
ANGSTROM_CUBED_TO_LITERS = 1 * 10 ** -27


class ConversionToolkit(Enum):
    """
    An enumeration of the different toolkits for converting MD files to other formats.
    """

    ParmEd = "parmed"
    InterMol = "intermol"
    # TopoTools = "topotools"


class PBCBox(Enum):
    """
    An enumeration of the different PBC box used when solvating with ``TLeap``.
    """
    cubic = "cubic"
    rectangular = "rectangular"
    octahedral = "octahedral"
