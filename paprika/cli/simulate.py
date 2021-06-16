import click


@click.command("simulate")
@click.option(
    "-s",
    "--script",
    "yaml_file",
    required=True,
    type=str,
    help="A YAML file containing the APR protocol for host-guest binding.",
)
def simulation():
    """Automatic pipeline for running APR calculations."""
