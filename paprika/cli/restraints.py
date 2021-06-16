import click

from .utils import OrderedGroup


@click.group()
def restraints():
    """Command-line interface for the `Restraints` modules."""
