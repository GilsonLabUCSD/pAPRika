import logging

import simtk.openmm as mm
import simtk.openmm.app as app
from simtk import unit

from parmed.openmm.reporters import NetCDFReporter

logger = logging.getLogger(__name__)


class OpenMM_GB_simulation:
    """Setup and run a GB simulation in OpenMM."""

    def __init__(self):

        self.path = "./"
        self.coordinates = None
        self.topology = None
        self.phase = None
        self.window = None

        self.min = dict(
            coordinates=self.coordinates,
            prefix="minimize",
            forcefield=None,
            platform="CUDA",
            devices="1",
            precision="mixed",
            max_iterations=5000,
            reporter_frequency=1000,
            tolerance=1,
            solvent=app.HCT,
            salt=0.1 * unit.mole / unit.liter,
            nonbonded_method=app.NoCutoff,
            nonbonded_cutoff=None,
            constraints=app.HBonds,
            temperature=300 * unit.kelvin,
            friction=1.0 / unit.picoseconds,
            timestep=2.0 * unit.femtoseconds,
            soft=False,
        )
        self.min["output"] = self.path + self.min["prefix"] + ".pdb"

        self.md = dict(
            coordinates=self.coordinates,
            minimized_coordinates=None,
            integrator=None,
            prefix="md",
            forcefield=None,
            platform="CUDA",
            devices="1",
            precision="mixed",
            reporter_frequency=1000,
            solvent=app.HCT,
            salt=0.1 * unit.mole / unit.liter,
            nonbonded_method=app.NoCutoff,
            nonbonded_cutoff=None,
            constraints=app.HBonds,
            temperature=300 * unit.kelvin,
            friction=1.0 / unit.picoseconds,
            timestep=2.0 * unit.femtoseconds,
            steps=10000,
        )
        self.md["output"] = self.path + self.md["prefix"] + ".nc"
        self.md["data"] = self.path + self.md["prefix"] + ".csv"

        def setup_system(self, settings, seed=None):
            """
            Provide a way to create an OpenMM system object with minimization or MD settings.
            Parameters
            ----------
            settings : dict
                A dictionary containing simulation settings.
            Returns
            -------
            simulation : simtk.openmm.app.Simulation
                The simulation object.
            system : simtk.openmm.System
                The system object.
            """

            prmtop = app.AmberPrmtopFile(self.topology)
            app.AmberInpcrdFile(settings["coordinates"])
            self.integrator = mm.LangevinIntegrator(
                settings["temperature"], settings["friction"], settings["timestep"]
            )

            if seed is not None:
                self.integrator.setRandomNumberSeed(seed)
            try:
                logger.debug(
                    "Integrator random number seed: {}".format(
                        self.integrator.getRandomNumberSeed()
                    )
                )
            except AttributeError:
                pass
            system = prmtop.createSystem(
                nonbondedMethod=settings["nonbonded_method"],
                implicitSolvent=settings["solvent"],
                implicitSolventSaltConc=settings["salt"],
                constraints=settings["constraints"],
            )

    return system