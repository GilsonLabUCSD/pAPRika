"""
This class is going to contain the simulation setup...
"""

import pkg_resources
import shutil
import textwrap
import parmed as pmd

from pathlib import Path
from paprika.tleap import System
from paprika.restraints import static_DAT_restraint, DAT_restraint
from paprika.restraints.restraints import create_window_list
from paprika.restraints.read_yaml import read_yaml
from paprika.io import save_restraints

import logging

logger = logging.getLogger(__name__)


def _get_installed_benchmarks():
    _installed_benchmarks = {}

    for entry_point in pkg_resources.iter_entry_points(group="taproom.benchmarks"):
        _installed_benchmarks[entry_point.name] = entry_point.load()
    return _installed_benchmarks


class Setup(object):
    """
    The Setup class provides a wrapper function around the preparation of the host-guest system and the application of restraints.
    """

    def __init__(self, host, guest, backend="openmm"):
        self.host = host
        self.guest = guest
        self.backend = backend

        self.directory = Path("benchmarks").joinpath(self.host).joinpath(self.guest)
        self.directory.mkdir(parents=True, exist_ok=True)

        if self.backend == "amber":
            raise NotImplementedError

        installed_benchmarks = get_benchmarks()
        host_yaml, guest_yaml = self.parse_yaml(installed_benchmarks)
        self.benchmark_path = host_yaml.parent
        self.host_yaml = read_yaml(host_yaml)
        self.guest_yaml = read_yaml(guest_yaml)

        self.build_bound()
        self.static_restraints,  self.conformational_restraints, self.wall_restraints,  self.guest_restraints = self.initialize_restraints()
        self.build_windows()

    def parse_yaml(self, installed_benchmarks):
        """
        Read the YAML recipe for the host and guest.

        Returns
        -------

        """
        try:
            host_yaml = installed_benchmarks["host_guest_pairs"][self.host]["yaml"]
        except KeyError:
            logger.error(f"Cannot find YAML recipe for host: {self.host}")
            logger.debug(installed_benchmarks)
            raise FileNotFoundError
        try:
            guest_yaml = installed_benchmarks["host_guest_pairs"][self.host][self.guest]
        except KeyError:
            logger.error(f"Cannot find YAML recipe for guest: {self.guest}")
            logger.debug(installed_benchmarks)
            raise FileNotFoundError
        return host_yaml, guest_yaml

    def build_bound(self):
        """
        This is going to change to accommodate SMIRNOFF99Frosst.

        Returns
        -------
        """
        system = System()
        system.output_path = self.directory
        system.output_prefix = f"{self.host}-{self.guest}"
        system.pbc_type = None
        system.neutralize = False
        system.template_lines = [
            "source leaprc.gaff",
            "source leaprc.water.tip3p",
            f"loadamberparams {self.benchmark_path.joinpath(self.host_yaml['frcmod'])}",
            f"loadamberparams {self.benchmark_path.joinpath('dummy.frcmod')}",
            f"{self.host.upper()} = loadmol2 {self.benchmark_path.joinpath(self.host_yaml['structure'])}",
            f"{self.guest.upper()} = loadmol2 {self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['structure'])}",
            f"DM1 = loadmol2 {self.benchmark_path.joinpath('dm1.mol2')}",
            f"DM2 = loadmol2 {self.benchmark_path.joinpath('dm2.mol2')}",
            f"DM3 = loadmol2 {self.benchmark_path.joinpath('dm3.mol2')}",
            f"model = loadpdb {self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['complex'])}",
        ]
        system.build()

    def initialize_restraints(self):

        windows = [
            self.host_yaml["calculation"]["windows"]["attach"],
            self.host_yaml["calculation"]["windows"]["pull"],
            None,
        ]
        structure = self.directory.joinpath(f"{self.host}-{self.guest}.pdb")

        static_restraints = []
        for restraint in self.host_yaml["calculation"]["restraints"]["static"]:
            static = static_DAT_restraint(
                restraint_mask_list=restraint["restraint"]["atoms"].split(),
                num_window_list=windows,
                ref_structure=str(structure),
                force_constant=restraint["restraint"]["force_constant"],
                amber_index=False if self.backend == "openmm" else True,
            )
            static_restraints.append(static)

        conformational_restraints = []
        if self.host_yaml["calculation"]["restraints"]["conformational"]:
            for conformational in self.host_yaml["calculation"]["restraints"][
                "conformational"
            ]:
                raise NotImplementedError
        else:
            logger.debug("Skipping conformational restraints...")

        wall_restraints = []
        if self.host_yaml["calculation"]["restraints"]["wall"]:
            for wall in self.host_yaml["calculation"]["restraints"]["wall"]:
                raise NotImplementedError
        else:
            logger.debug("Skipping wall restraints...")

        guest_restraints = []
        for restraint in self.guest_yaml["restraints"]:
            mask = restraint["restraint"]["atoms"].split()

            guest_restraint = DAT_restraint()
            guest_restraint.auto_apr = True
            guest_restraint.continuous_apr = True
            guest_restraint.amber_index = False if self.backend == "openmm" else True
            guest_restraint.topology = str(structure)
            guest_restraint.mask1 = mask[0]
            guest_restraint.mask2 = mask[1]
            guest_restraint.mask3 = mask[2] if len(mask) > 2 else None
            guest_restraint.mask4 = mask[3] if len(mask) > 3 else None

            guest_restraint.attach["target"] = restraint["restraint"]["attach"][
                "target"
            ]
            guest_restraint.attach["fc_final"] = restraint["restraint"]["attach"][
                "force_constant"
            ]
            guest_restraint.attach["fraction_list"] = self.host_yaml["calculation"][
                "lambda"
            ]["attach"]

            guest_restraint.pull["target_final"] = self.host_yaml["calculation"][
                "target"
            ]["pull"]
            guest_restraint.pull["num_windows"] = windows[1]

            guest_restraint.initialize()
            guest_restraints.append(guest_restraint)

        return static_restraints,        conformational_restraints,        wall_restraints,         guest_restraints

    def build_windows(self):
        save_restraints(self.static_restraints + self.conformational_restraints + self.wall_restraints + self.guest_restraints, self.directory.joinpath("restraints.json"))
        window_list = create_window_list(self.guest_restraints)
        
        for window in window_list:
            self.directory.joinpath("windows").joinpath(window).mkdir(parents=True, exist_ok=True)
            self.translate(window)
            self.solvate(window)

            if self.backend == "amber":
                # Write the restraint file in each window.
                raise NotImplementedError
            else:
                self.initialize_calculation(window)
    
    def translate(self, window):
        window_path = self.directory.joinpath("windows").joinpath(window)
        if window[0] == "a":
            # Copy the initial structure.
            source_prmtop = self.directory.joinpath(f"{self.host}-{self.guest}.prmtop")
            source_inpcrd = self.directory.joinpath(f"{self.host}-{self.guest}.rst7")
            shutil.copy(source_prmtop, window_path)
            shutil.copy(source_inpcrd, window_path)
        elif window[0] == "p":
            # Translate the guest.
            source_prmtop = self.directory.joinpath(f"{self.host}-{self.guest}.prmtop")
            source_inpcrd = self.directory.joinpath(f"{self.host}-{self.guest}.rst7")

            structure = pmd.load_file(str(source_prmtop), str(source_inpcrd), structure = True)
            target_difference = self.guest_restraints[0].phase['pull']['targets'][int(window[1:])] - self.guest_restraints[0].pull['target_initial']
            for atom in structure.atoms:
                if atom.residue.name == self.guest.upper():
                    atom.xz += target_difference
            structure.save(str(window_path.joinpath(f"{self.host}-{self.guest}.prmtop")), overwrite=True)
            structure.save(str(window_path.joinpath(f"{self.host}-{self.guest}.rst7")))

        elif window[0] == "r":
            # Copy the final pull window.
            source_prmtop = self.directory.joinpath("windows").joinpath(f"p{self.host_yaml['calculation']['windows']['pull']:03d}").joinpath(f"{self.host}-{self.guest}.prmtop")
            source_inpcrd = self.directory.joinpath("windows").joinpath(f"p{self.host_yaml['calculation']['windows']['pull']:03d}").joinpath(f"{self.host}-{self.guest}.rst7")
            shutil.copy(source_prmtop, window_path)
            shutil.copy(source_inpcrd, window_path)

    def solvate(self, window):
        window_path = self.directory.joinpath("windows").joinpath(window)
        source_prmtop = window_path.joinpath(f"{self.host}-{self.guest}.prmtop")
        source_inpcrd = window_path.joinpath(f"{self.host}-{self.guest}.rst7")
        pdb = window_path.joinpath(f"{self.host}-{self.guest}.pdb")
        
        if not Path.exists(pdb):
            structure = pmd.load_file(str(source_prmtop), str(source_inpcrd), structure = True)
            structure.save(str(pdb), overwrite=True)
    
        system = System()
        system.output_path = window_path
        system.output_prefix = f"{self.host}-{self.guest}-sol"
        
        system.target_waters = self.host_yaml["calculation"]["waters"]
        system.neutralize = True
        # system.add_ions = ["Na+", 6, "Cl-", 6]
        system.template_lines = [
        "source leaprc.gaff",
        "source leaprc.water.tip3p",
        f"loadamberparams {self.benchmark_path.joinpath(self.host_yaml['frcmod'])}",
        f"loadamberparams {self.benchmark_path.joinpath('dummy.frcmod')}",
        f"{self.host.upper()} = loadmol2 {self.benchmark_path.joinpath(self.host_yaml['structure'])}",
        f"{self.guest.upper()} = loadmol2 {self.benchmark_path.joinpath(self.guest).joinpath(self.guest_yaml['structure'])}",
        f"DM1 = loadmol2 {self.benchmark_path.joinpath('dm1.mol2')}",
        f"DM2 = loadmol2 {self.benchmark_path.joinpath('dm2.mol2')}",
        f"DM3 = loadmol2 {self.benchmark_path.joinpath('dm3.mol2')}",
        f"model = loadpdb {self.host}-{self.guest}.pdb",
        ]
        system.build()

    def initialize_calculation(self, window):
        if self.backend == "amber":
            # Write simulation input files in each directory
            raise NotImplementedError

        window_path = self.directory.joinpath("windows").joinpath(window)
        source_prmtop = window_path.joinpath(f"{self.host}-{self.guest}.prmtop")
        source_inpcrd = window_path.joinpath(f"{self.host}-{self.guest}.rst7")

        
        openmm_string = f"""
        import os
        import re
        import shutil
        import json

        import numpy as np
        import parmed as pmd

        import simtk.openmm as mm
        import simtk.openmm.app as app
        import simtk.unit as unit
        from parmed.openmm.reporters import NetCDFReporter

        settings = {{
        'nonbonded_method': app.PME,
        'nonbonded_cutoff': 8.0 * unit.angstrom,
        'temperature': 300 * unit.kelvin,
        'timestep': 0.002 * unit.picosecond,
        'constraints': app.HBonds,
        'friction': 1.0 / unit.picoseconds,
        }}

        prmtop_file = "{self.host}-{self.guest}-sol.prmtop"
        prmtop = app.AmberPrmtopFile(prmtop_file)

        system = prmtop.createSystem(
                nonbondedMethod=settings["nonbonded_method"],
                nonbondedCutoff=settings["nonbonded_cutoff"],
                constraints=settings["constraints"],
        )

        integrator = mm.LangevinIntegrator(
                    settings["temperature"],
                    settings["friction"],
                    settings["timestep"]
        )
        """
        for restraint in self.static_restraints + self.conformational_restraints + self.guest_restraints:
            openmm_string += openmm_string_from_restraint(restraint, window)
        for restraint in self.wall_restraints:
            # Custom force here...
            raise NotImplementedError

        openmm_string += f"""
        system.addForce(mm.MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin, 100))

        inpcrd = app.AmberInpcrdFile("{self.host}-{self.guest}-sol.rst7")
        simulation = app.Simulation(prmtop.topology, system, integrator,
                                mm.Platform.getPlatformByName('CUDA'))

        simulation.context.setPositions(inpcrd.positions)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(settings["temperature"] * unit.kelvin)

        simulation.reporters.append(NetCDFReporter(f"openmm.nc", 500))
        simulation.reporters.append(app.StateDataReporter("openmm.log", 5000, step=True,
            potentialEnergy=True, temperature=True, speed=True))
        simulation.step(5000000)
        """
        with open(window_path.joinpath("openmm.py"), "w") as f:
            f.write(textwrap.dedent(openmm_string))

def get_benchmarks():
    """
    Determine the installed benchmarks.

    """
    installed_benchmarks = _get_installed_benchmarks()
    return installed_benchmarks

def openmm_string_from_restraint(restraint, window):
    if window[0] == "a":
        phase = "attach"
    elif window[0] == "p":
        phase = "pull"
    elif window[0] == "r":
        phase = "release"
    window_number = int(window[1:])

    if restraint.mask2 and not restraint.mask3:
        if not restraint.group1 and not restraint.group2:
            string = f"""
        bond_restraint = mm.CustomBondForce('k * (r - r_0)^2')
        bond_restraint.addPerBondParameter('k')
        bond_restraint.addPerBondParameter('r_0')

        r_0 = {restraint.phase[phase]["targets"][window_number]} * unit.angstroms
        k = {restraint.phase[phase]["force_constants"][window_number]} * unit.kilocalories_per_mole / unit.angstroms**2 
        bond_restraint.addBond({restraint.index1[0]}, {restraint.index2[0]}, 
        [k, r_0])
        system.addForce(bond_restraint)
            """
        else:
            # Probably needs mm.CustomCentroidBondForce (?)
            raise NotImplementedError
    elif restraint.mask3 and not restraint.mask4:
        if not restraint.group1 and not restraint.group2 and not restraint.group3:
            string = f"""
        angle_restraint = mm.CustomAngleForce('k * (theta - theta_0)^2')
        angle_restraint.addPerAngleParameter('k')
        angle_restraint.addPerAngleParameter('theta_0')

        theta_0 = {restraint.phase[phase]["targets"][window_number]} * unit.degrees
        k = {restraint.phase[phase]["force_constants"][window_number]} * unit.kilocalories_per_mole / unit.degrees**2 
        angle_restraint.addAngle({restraint.index1[0]}, {restraint.index2[0]}, {restraint.index3[0]}, 
        [k, theta_0])
        system.addForce(angle_restraint)
            """
        else:
            # Probably needs mm.CustomCentroidAngleForce (?)
            raise NotImplementedError

    elif restraint.mask4:
        if not restraint.group1 and not restraint.group2 and not restraint.group3 and not restraint.group4:
            string = f"""
        dihedral_restraint = mm.CustomTorsionForce('k * (theta - theta_0)^2')
        dihedral_restraint.addPerTorsionParameter('k')
        dihedral_restraint.addPerTorsionParameter('theta_0')

        theta_0 = {restraint.phase[phase]["targets"][window_number]} * unit.degrees
        k = {restraint.phase[phase]["force_constants"][window_number]} * unit.kilocalories_per_mole / unit.degrees**2 
        dihedral_restraint.addTorsion({restraint.index1[0]}, {restraint.index2[0]}, {restraint.index3[0]}, {restraint.index4[0]},
        [k, theta_0])
        system.addForce(dihedral_restraint)
            """
        else:
            # Probably needs mm.CustomCentroidTorsionForce (?)
            raise NotImplementedError

    return string
