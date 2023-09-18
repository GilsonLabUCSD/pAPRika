import os
import subprocess
from typing import List

import openmm.unit as openmm_unit
import parmed
from openff.toolkit.topology import Molecule
from openff.units import unit
from openff.units.openmm import to_openmm
from openmm.app import PDBFile
from rdkit.Chem import AllChem as Chem


class CyclodextrinFragmenter:
    """
    The class that defines a Cyclodextrin molecule with methods to assign partial charges
    using a fragment-based approach. This is adapted from Tobias Hüfner's code:

    https://github.com/wutobias/collection/blob/master/python/am1bcc-glyco-fragmenter.py
    """

    @property
    def monomers(self) -> List[Chem.Mol]:
        """list: A list of monomers."""
        return self._monomers

    @monomers.setter
    def monomers(self, value: List[Chem.Mol]):
        self._monomers = value

    @property
    def input_molecule(self) -> Chem.Mol:
        """Chem.Mol: The Cyclodextrin molecule."""
        return self._input_molecule

    @input_molecule.setter
    def input_molecule(self, value: Chem.Mol):
        self._input_molecule = value

    @property
    def dummy_atom_labels(self) -> List[int]:
        """list: List of dummy atom labels -- for O1 and C4 bridge."""
        return self._dummy_atom_labels

    @property
    def label_O1(self) -> int:
        """int: Integer label for atom O1"""
        return self._label_O1

    @label_O1.setter
    def label_O1(self, value: int):
        self._label_O1 = value
        self._dummy_atom_labels[0] = value

    @property
    def label_C4(self) -> int:
        """int: Integer label for atom C4"""
        return self._label_C4

    @label_C4.setter
    def label_C4(self, value: int):
        self._label_C4 = value
        self._dummy_atom_labels[1] = value

    def __init__(self, input_molecule: Chem.Mol, monomers: List[Chem.Mol] = None):
        self._input_molecule = input_molecule
        self._monomers = [] if monomers is None else monomers
        self._capped_monomers = []
        self._capped_monomers_atom_indices = []
        self._capped_monomers_smiles = []
        self._capped_monomers_index = []
        self._capped_monomers_charges = {}
        self._label_O1 = 1
        self._label_C4 = 4
        self._dummy_atom_labels = [self._label_O1, self._label_C4]
        self._partial_charge_list = []
        self._partial_charge_dictionary = {}

    def add_monomer(self, molecule: Chem.Mol):
        """Utility function that adds monomers to the list."""
        self.monomers.append(molecule)

    def assign_partial_charges(self, partial_charge_method="am1bcc"):
        """
        Assign partial charges to cyclodextrin molecule using a fragment-based approach

        Parameters
        ----------
        partial_charge_method: str
            The partial charge method (see OpenFF-Toolkit documentation)

        """

        # Generate charges for each capped monomer
        openff_list = dict()
        self._capped_monomers_charges = dict()
        self._get_capped_monomer()

        for fragment_index in list(set(self._capped_monomers_index)):
            # Generate Conformers
            molecule = Chem.AddHs(self._capped_monomers[fragment_index])
            Chem.EmbedMultipleConfs(molecule)
            Chem.MMFFOptimizeMoleculeConfs(molecule)

            # Write SDF file for methyl-capped fragment
            writer = Chem.SDWriter(f"./frag{fragment_index}.sdf")
            writer.write(molecule, confId=1)
            writer.close()

            # Assign partial charges using OFF-Toolkit with Antechamber Backend
            molecule_openff = Molecule.from_file(f"./frag{fragment_index}.sdf")
            molecule_openff.assign_partial_charges(
                partial_charge_method=partial_charge_method,
                normalize_partial_charges=True,
            )

            # Add molecule and partial charges to list
            openff_list[fragment_index] = molecule_openff
            try:
                self._capped_monomers_charges[
                    fragment_index
                ] = molecule_openff.partial_charges.m_as(unit.elementary_charge)
            except AttributeError:
                self._capped_monomers_charges[
                    fragment_index
                ] = molecule_openff.partial_charges.value_in_unit(
                    openmm_unit.elementary_charge
                )

            # Clean up
            os.remove(f"./frag{fragment_index}.sdf")

        # Add monomer charges to cyclodextrin molecule
        self._partial_charge_dictionary = dict()
        current_charge = 0.0

        # Loop over monomers
        for monomer_index, fragment_index in enumerate(self._capped_monomers_index):
            molecule_atom_indices = self._capped_monomers_atom_indices[monomer_index]
            molecule = self._capped_monomers[monomer_index]
            self_match = (
                openff_list[fragment_index].to_rdkit().GetSubstructMatch(molecule)
            )
            for (
                fragment_atom_index,
                molecule_atom_index,
            ) in enumerate(molecule_atom_indices):
                partial_charge = self._capped_monomers_charges[fragment_index][
                    self_match[fragment_atom_index]
                ]
                self._partial_charge_dictionary[molecule_atom_index] = partial_charge
                current_charge += partial_charge

        # Normalize partial charges
        expected_charge = float(Chem.GetFormalCharge(self.input_molecule))
        charge_offset = (current_charge - expected_charge) / float(
            len(self._partial_charge_dictionary)
        )
        self._partial_charge_list = list()

        for atom_index in range(self.input_molecule.GetNumAtoms()):
            # Skip capped atoms
            if atom_index not in self._partial_charge_dictionary:
                continue

            charge = self._partial_charge_dictionary[atom_index] - charge_offset
            self._partial_charge_list.append(str(charge))

        # Add partial charges to RDKit Molecule instance
        self.input_molecule.SetProp(
            "atom.dprop.PartialCharge", "\n".join(self._partial_charge_list)
        )

    def to_file(self, file, file_format, residue_name="MGO", atom_type="sybyl"):
        """
        Save the cyclodextrin molecule to file.

        Parameters
        ----------
        file: str
            File name
        file_format: str
            File format -- either SDF, MOL2 and PDB (through parmed)
        residue_name: str
            Residue name -- mainly for MOL2 files
        atom_type: str
            The atom type in MOL2 files -- gaff, gaff2, amber, bcc or sybyl
        """
        off_molecule = Molecule.from_rdkit(self.input_molecule)

        # Set residue name
        for atom in off_molecule.atoms:
            atom.metadata["residue_name"] = residue_name

        # Save to file
        if file_format == "SDF":
            off_molecule.to_file(file, file_format=file_format)

        elif file_format == "PDB":
            with open(file, "w") as f:
                PDBFile.writeFile(
                    off_molecule.to_topology().to_openmm(),
                    to_openmm(off_molecule.conformers[0]),
                    f,
                    keepIds=True,
                )

        elif file_format == "MOL2":
            # Convert SDF to MOL2 as boiler plate
            off_molecule.to_file("temp.sdf", file_format="SDF")
            output = subprocess.Popen(
                [
                    "antechamber",
                    "-fi",
                    "sdf",
                    "-fo",
                    "mol2",
                    "-i",
                    "temp.sdf",
                    "-o",
                    "temp1.mol2",
                    "-rn",
                    residue_name,
                    "-at",
                    "sybyl",
                    "-pf",
                    "y",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd="./",
            )
            output = output.stdout.read().decode().splitlines()

            # Load temp.mol2 with SYBYL atom types and add charges to atoms
            structure = parmed.load_file("temp1.mol2", structure=True)
            for off_atom, atom, charge in zip(
                off_molecule.atoms, structure.atoms, self._partial_charge_list
            ):
                atom.charge = float(charge)
                atom.name = off_atom.name

            # Save current structure to MOL2 file
            structure.save("temp2.mol2", overwrite=True)

            # Rewrite with residue name and correct atom type using Antechamber
            output = subprocess.Popen(
                [
                    "antechamber",
                    "-fi",
                    "mol2",
                    "-fo",
                    "mol2",
                    "-i",
                    "temp2.mol2",
                    "-o",
                    file,
                    "-rn",
                    residue_name,
                    "-at",
                    atom_type.lower(),
                    "-pf",
                    "y",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd="./",
            )
            output = output.stdout.read().decode().splitlines()

            # Clean-up temp file
            os.remove("./temp.sdf")
            os.remove("./temp1.mol2")
            os.remove("./temp2.mol2")

    @staticmethod
    def _get_atom_map_number(atom: Chem.Atom) -> int:
        """
        Given a RDKit Atom, return the atom map number for an atom.
        For cyclodextrin, I labelled `1` and `4` for O1 and C4 atoms, respectively.
        """
        if "molAtomMapNumber" in atom.GetPropNames():
            i = int(atom.GetProp("molAtomMapNumber"))
        else:
            i = -1

        return i

    def _get_capped_monomer(self):
        """
        Get matched monomers and cap terminals with methyl.
        """
        self._capped_monomers = list()
        self._capped_monomers_atom_indices = list()

        # Get bond list
        bond_list, dummy_labels = self._get_bond_list()

        # Bond break based on bond list and dummy labels
        input_molecule_fragment = Chem.FragmentOnBonds(
            self.input_molecule,
            bond_list,
            dummyLabels=dummy_labels,
        )

        # Loop over fragments
        for fragment_molecule, fragment_atom_indices in zip(
            Chem.GetMolFrags(input_molecule_fragment, asMols=True),
            Chem.GetMolFrags(input_molecule_fragment, asMols=False),
        ):
            atom_index_list = list()
            rwmol = Chem.RWMol(fragment_molecule)

            for atom, atom_index in zip(
                fragment_molecule.GetAtoms(), fragment_atom_indices
            ):
                # Get label
                label = atom.GetIsotope()

                # Add Methyl to O1
                if label == self.label_O1:
                    rwmol.ReplaceAtom(atom.GetIdx(), Chem.Atom(8))
                    rwmol.AddAtom(Chem.Atom(6))
                    rwmol.AddBond(
                        atom.GetIdx(),
                        fragment_molecule.GetNumAtoms(),
                        order=Chem.BondType.SINGLE,
                    )
                    fragment_molecule = rwmol.GetMol()

                # Add Methyl to C4
                elif label == self.label_C4:
                    rwmol.ReplaceAtom(atom.GetIdx(), Chem.Atom(6))

                # Add atom index to list
                else:
                    atom_index_list.append(atom_index)

            # Generate 2D coordinates
            fragment_molecule = rwmol.GetMol()
            Chem.SanitizeMol(fragment_molecule)
            fragment_molecule.Compute2DCoords()

            # Add capped fragment molecule to list
            self._capped_monomers.append(fragment_molecule)
            self._capped_monomers_atom_indices.append(atom_index_list)

        # Generate unique smiles and index
        smiles_list = [
            Chem.CanonSmiles(Chem.MolToSmiles(mol)) for mol in self._capped_monomers
        ]
        self._capped_monomers_smiles = list(set(smiles_list))
        self._capped_monomers_index = [
            smiles_list.index(smiles) for smiles in smiles_list
        ]

    def _get_bond_list(self):
        """
        Get bond index between monomers - used to cut the macrocylic molecule
        """
        bond_list = list()
        dummy_labels = list()

        # Loop over monomers
        for monomer in self.monomers:
            # Find monomer in molecule
            input_matches = self.input_molecule.GetSubstructMatches(monomer)

            # Loop over monomer matches
            for match in input_matches:
                # Loop over atom indices in `match`
                for query_index, atom_index in enumerate(match):
                    # Select only O1 atom.
                    if (
                        self._get_atom_map_number(monomer.GetAtomWithIdx(query_index))
                        != 1
                    ):
                        continue

                    atom = self.input_molecule.GetAtomWithIdx(atom_index)

                    # Loop over neighbors
                    for neighbor_atom in atom.GetNeighbors():
                        # Select neighboring atom not in `SubstructMatch` (or in monomer)
                        if neighbor_atom.GetIdx() in match:
                            continue

                        # Get bond between atom O1 and its neighbor (i.e. C4)
                        bond = self.input_molecule.GetBondBetweenAtoms(
                            atom_index, neighbor_atom.GetIdx()
                        )
                        bond_index = bond.GetIdx()

                        # Add to list unique bond_index
                        if bond_index in bond_list:
                            continue

                        bond_list.append(bond_index)

                        # Set labels for dummy
                        if bond.GetBeginAtomIdx() == atom_index:
                            dummy_labels.append(self.dummy_atom_labels)
                        else:
                            dummy_labels.append(self.dummy_atom_labels[::-1])

        return bond_list, dummy_labels
