"""
Functions for generate GAFF files and parameters.
"""

import logging
import os
import subprocess
from typing import Optional

import numpy as np
import parmed

logger = logging.getLogger(__name__)
_PI_ = np.pi


def generate_gaff(
    mol2_file: str,
    residue_name: str,
    output_name: Optional[str] = None,
    need_gaff_atom_types: Optional[bool] = True,
    generate_frcmod: Optional[bool] = True,
    directory_path: Optional[str] = "benchmarks",
    gaff_version: Optional[str] = "gaff2",
):
    """
    Module to generate GAFF files given a mol2 file.

    Parameters
    ----------
    mol2_file
        The name of the mol2 structure file.
    residue_name
        The residue name of the molecule.
    output_name
        The name for the output file.
    need_gaff_atom_types
        Whether to generate GAFF atoms or not. Currently, this is the only choice.
    gaff_version
        The GAFF version to use ("gaff1", "gaff2")
    generate_frcmod
        Option to generate a GAFF frcmod file.
    directory_path
        The working directory where the files will be stored.
    """

    if output_name is None:
        output_name = mol2_file.stem

    if need_gaff_atom_types:
        _generate_gaff_atom_types(
            mol2_file=mol2_file,
            residue_name=residue_name,
            output_name=output_name,
            gaff_version=gaff_version,
            directory_path=directory_path,
        )
        logging.debug(
            "Checking to see if we have a multi-residue MOL2 file that should be converted "
            "to single-residue..."
        )
        structure = parmed.load_file(
            os.path.join(directory_path, f"{output_name}.{gaff_version}.mol2"),
            structure=True,
        )
        if len(structure.residues) > 1:
            structure[":1"].save("tmp.mol2")
            if os.path.exists("tmp.mol2"):
                os.rename(
                    "tmp.mol2",
                    os.path.join(directory_path, f"{output_name}.{gaff_version}.mol2"),
                )
                logging.debug("Saved single-residue MOL2 file for `tleap`.")
            else:
                raise RuntimeError(
                    "Unable to convert multi-residue MOL2 file to single-residue for `tleap`."
                )

        if generate_frcmod:
            _generate_frcmod(
                mol2_file=f"{output_name}.{gaff_version}.mol2",
                gaff_version=gaff_version,
                output_name=output_name,
                directory_path=directory_path,
            )
        else:
            raise NotImplementedError()

    else:
        raise NotImplementedError()


def _generate_gaff_atom_types(
    mol2_file: str,
    residue_name: str,
    output_name: str,
    gaff_version: Optional[str] = "gaff2",
    directory_path: Optional[str] = "benchmarks",
):
    """Generate a mol2 file with GAFF atom types."""

    if gaff_version.lower() not in ["gaff", "gaff2"]:
        raise KeyError(
            f"Parameter set {gaff_version} not supported. Only [gaff, gaff2] are allowed."
        )

    p = subprocess.Popen(
        [
            "antechamber",
            "-i",
            mol2_file,
            "-fi",
            "mol2",
            "-o",
            f"{output_name}.{gaff_version}.mol2",
            "-fo",
            "mol2",
            "-rn",
            f"{residue_name.upper()}",
            "-at",
            f"{gaff_version}",
            "-an",
            "no",
            "-dr",
            "no",
            "-pf",
            "yes",
        ],
        cwd=directory_path,
    )
    p.communicate()
    print(p)

    remove_files = [
        "ANTECHAMBER_AC.AC",
        "ANTECHAMBER_AC.AC0",
        "ANTECHAMBER_BOND_TYPE.AC",
        "ANTECHAMBER_BOND_TYPE.AC0",
        "ATOMTYPE.INF",
    ]
    files = [os.path.join(directory_path, file) for file in remove_files]
    for file in files:
        if os.path.isfile(file):
            logger.debug(f"Removing temporary file: {file}")
            file.unlink()

    if not os.path.exists(f"{output_name}.{gaff_version}.mol2"):
        # Try with the newer (AmberTools 19) version of `antechamber` which doesn't have the `-dr` flag
        p = subprocess.Popen(
            [
                "antechamber",
                "-i",
                mol2_file,
                "-fi",
                "mol2",
                "-o",
                f"{output_name}.{gaff_version}.mol2",
                "-fo",
                "mol2",
                "-rn",
                f"{residue_name.upper()}",
                "-at",
                f"{gaff_version}",
                "-an",
                "no",
                "-pf",
                "yes",
            ],
            cwd=directory_path,
        )
        p.communicate()

        remove_files = [
            "ANTECHAMBER_AC.AC",
            "ANTECHAMBER_AC.AC0",
            "ANTECHAMBER_BOND_TYPE.AC",
            "ANTECHAMBER_BOND_TYPE.AC0",
            "ATOMTYPE.INF",
        ]
        files = [os.path.join(directory_path, file) for file in remove_files]
        for file in files:
            if os.path.isfile(file):
                logger.debug(f"Removing temporary file: {file}")
                file.unlink()


def _generate_frcmod(
    mol2_file: str,
    gaff_version: str,
    output_name: str,
    directory_path: Optional[str] = "benchmarks",
):
    """Generate an AMBER .frcmod file given a mol2 file."""

    if gaff_version.lower() not in ["gaff", "gaff2"]:
        raise KeyError(
            f"Parameter set {gaff_version} not supported. Only [gaff, gaff2] are allowed."
        )

    subprocess.Popen(
        [
            "parmchk2",
            "-i",
            str(mol2_file),
            "-f",
            "mol2",
            "-o",
            f"{output_name}.{gaff_version}.frcmod",
            "-s",
            f"{gaff_version}",
        ],
        cwd=directory_path,
    )
