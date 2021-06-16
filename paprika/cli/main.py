import click

from .analysis import analyze
from .build import align, dummy_atoms
from .restraints import restraints
from .simulate import simulation
from .utils import OrderedGroup


@click.group(cls=OrderedGroup)
def main():
    """pAPRika command-line interface."""


# Add cli groups
main.add_command(analyze)
main.add_command(simulation)
main.add_command(align)
main.add_command(dummy_atoms)
main.add_command(restraints)


if __name__ == "__main__":
    main()
