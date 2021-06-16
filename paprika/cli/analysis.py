import click

from .utils import OrderedGroup


@click.group(cls=OrderedGroup)
@click.option(
    "-top",
    "--topology",
    "topology",
    required=True,
    type=str,
    help="The topology used in the simulation.",
)
@click.option(
    "-traj",
    "--trajectory",
    "trajectory",
    required=True,
    type=str,
    help="The trajectory files, can include wildcards.",
)
@click.option(
    "-restr",
    "--restraints",
    "restraints_file",
    required=True,
    type=str,
    help="The JSON file containing the restraints definitions.",
)
@click.option(
    "-nboot",
    "--boot_cycles",
    "boot_cycles",
    default=1000,
    type=int,
    help="The number of iterations for bootstrapping.",
)
@click.option(
    "--output",
    "output_file",
    default="results.json",
    type=str,
    help="The name of the output file (default: results.json).",
)
@click.option(
    "-T",
    "--temperature",
    "temperature",
    default=298.15,
    type=float,
    help="The temperature of the system (default: 298.15 Kelvin).",
)
@click.option(
    "--ti_matrix",
    "ti_matrix",
    default="diagonal",
    type=str,
    help="Option to estimate the mean and SEM in TI: 'full', 'diagonal', 'endpoints' (default: diagonal)",
)
@click.option(
    "-roi",
    "--compute_roi",
    "compute_roi",
    is_flag=True,
    help="Option to compute the return on investment (ROI) value for each window for TI",
)
@click.option(
    "--methods",
    "methods",
    nargs=0,
    default=["ti-block"],
    type=str,
    help="The method to estimate the free energy ('ti-block', 'mbar-block' or 'mbar-autoc') (default: ti-block)",
)
@click.argument("methods", nargs=-1)
def analyze(trajectory, boot_cycles):
    """Analyze MD trajectories to extract free energies."""
