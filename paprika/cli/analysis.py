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
    "-o",
    "--output",
    "output_file",
    default="results.json",
    type=str,
    help="The name of the output file (default: results.json).",
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
@click.option(
    "-N",
    "--boot_cycles",
    "boot_cycles",
    default=1000,
    type=int,
    help="The number of iterations for bootstrapping (default: 1000).",
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
    "--compute_roi",
    "compute_roi",
    is_flag=True,
    help="Option to compute the return on investment (ROI) value for each window in TI",
)
@click.option(
    "--conservative",
    "conservative_subsample",
    is_flag=True,
    help="Option to round the statistical inefficiency g to the nearest integer for MBAR."
)
def analyze(
    topology,
    trajectory,
    restraints_file,
    path,
    output_files,
    methods,
    boot_cycles,
    temperature,
    ti_matrix,
    compute_roi,
):
    """Analyze MD trajectories to extract free energies."""
    from paprika.analysis import fe_calc
    from paprika.io import load_restraints

    restraints = load_restraints(restraints_file)

    free_energy = fe_calc()
    free_energy.temperature = temperature
    free_energy.topology = topology
    free_energy.trajectory = trajectory
    free_energy.restraint_list = restraints
    free_energy.path = path

    free_energy.collect_data()

    free_energy.methods = methods
    free_energy.ti_matrix = ti_matrix
    free_energy.bootcycles = boot_cycles
    free_energy.compute_roi = compute_roi

    free_energy.compute_free_energy()
