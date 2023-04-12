import numpy as np
import pyunitwizard as puw
from rdkit import Chem

from typing import List

from openpharmacophore.molecular_systems.topology import Topology, CHAIN_NAMES, SOLVENT_AND_IONS
from openpharmacophore.molecular_systems.ligand import ligand_from_topology
from openpharmacophore.utils.maths import delete
from openpharmacophore.molecular_systems import convert


class Protein:
    """ Represents a protein and optionally its different structures
        if it comes from a MD trajectory. The protein may contain ligands
        but is not necessary.

        Parameters
        ----------
        topology : Topology
            The topology.

        coords : Quantity
            A quantity of shape (n_structures, n_atoms, 3)
    """
    def __init__(self, topology, coords):
        self._validate_coords(topology, coords)
        self._topology = topology
        self._coords = coords

    @property
    def n_atoms(self) -> int:
        return self._topology.n_atoms

    @property
    def n_residues(self) -> int:
        return self._topology.n_residues

    @property
    def n_chains(self) -> int:
        return self.topology.n_chains

    @property
    def topology(self) -> Topology:
        return self._topology

    @property
    def coords(self) -> np.ndarray:
        return self._coords

    def has_hydrogens(self) -> bool:
        return self._topology.has_hydrogens()

    def has_ligands(self) -> bool:
        return self._topology.has_ligands()

    def ligand_ids(self) -> List[str]:
        return self._topology.ligand_ids()

    def has_solvent_or_ions(self) -> bool:
        return self._topology.has_solvent_or_ions()

    def get_ligand(self, ligand_id, remove_hyd=False):
        """ Extract a ligand assuming there is one.

            Parameters
            ----------
            ligand_id : str
                ID of the ligand. ID is composed by pdb id followed by chain name,
                i.e. "EST:B".

            remove_hyd : bool
                Whether to remove hydrogens from the ligand.

            Returns
            -------
            openpharmacophore.Ligand
                A ligand. Its bond orders may be incorrect.
        """
        lig_name, chain = ligand_id.split(":")
        lig_indices = self._topology.get_residue_indices(lig_name, chain)
        assert len(lig_indices) > 0, f"Ligand {ligand_id} could not be extracted"

        lig_coords = self._coords[:, lig_indices, :]
        ligand_top = self._topology.subset(lig_indices)
        return ligand_from_topology(ligand_top, lig_coords, remove_hyd)

    def remove_ligand(self, ligand_id):
        """ Remove a ligand from this protein.

            Parameters
            ----------
            ligand_id : str
                ID of the ligand

        """
        lig_name, chain = ligand_id.split(":")
        lig_indices = self._topology.get_residue_indices(lig_name, chain)
        self._topology.remove_atoms(lig_indices, inplace=True)
        self._coords = delete(self._coords, lig_indices, axis=1)

    def remove_all_ligands(self):
        """ Remove all ligands from this protein.
        """
        for lig in self.ligand_ids():
            self.remove_ligand(lig)

    def add_hydrogens(self):
        """ Add hydrogens to the protein.

            Assumes that the protein has a single frame.
        """
        mol = convert.topology_to_mol(self.topology, puw.get_value(self.coords[0], "nanometers"))
        mol = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
        # Note: Hs added with rdkit residue number does not correspond to the original residue number

        self._topology = convert.mol_to_topology(mol)

        coords = puw.quantity(mol.GetConformer(0).GetPositions(), "angstroms")
        self._coords = np.expand_dims(coords, axis=0)

        self._validate_coords(self.topology, self._coords)

    def atoms_at_distance(self, frame, centroid, max_dist, min_dist=0):
        """ Get the indices of the atoms that are at
            min_dist <= centroid <= max_dist.

            Parameters
            ----------
            frame : int
                Frame of the trajectory
            centroid : puw.Quantity
                Shape (1,3)
            max_dist : puw.Quantity
                Scalar
            min_dist : puw.Quantity
                Scalar

            Returns
            -------
            np.ndarray
                Array if integers with the indices.

        """
        if min_dist == 0:
            min_dist = puw.quantity(0, "nanometers")

        distance = np.sqrt(np.sum(np.power(self._coords[frame] - centroid, 2), axis=1))
        return np.where(
            np.logical_and(distance > min_dist, distance <= max_dist)
        )[0]

    def slice(self, atoms, frame=None):
        """ Obtain a subset of the protein with the specified atoms.

            Parameters
            ----------
            atoms : list[int]
            atoms : list[int]
                Indices of the atoms in the sliced protein.

            frame : int, optional
                Frame to extract. If None all frames are extracted.

            Returns
            -------
            Protein
                The sliced protein.
        """
        topology = self._topology.subset(atoms)
        if frame is None:
            coords = self._coords[:, atoms, :]
        else:
            coords = self._coords[frame, atoms, :]
            coords = np.expand_dims(coords, axis=0)

        return Protein(topology, coords)

    def concatenate(self, topology, coords):
        """ Concatenate new residues to this protein.

            Parameters
            ----------
            topology : Topology

            coords : QuantityLike

        """
        self._topology = self._topology.join(topology)
        self._coords = np.concatenate((self._coords, coords), axis=1)
        self._validate_coords(self._topology, self._coords)

    def extract_chain(self, chain_id):
        """ Extract the chain with given id from the topology.

            Parameters
            ----------
            chain_id : str or int
                Id of the chain
        """
        if isinstance(chain_id, str):
            chain_id = CHAIN_NAMES.index(chain_id)

        atoms = []
        for at in self.topology.iter_atoms():
            if at.residue.chain.index == chain_id:
                atoms.append(at.index)

        self._topology = self._topology.subset(atoms)
        self._coords = self.coords[:, atoms, :]
        self._validate_coords(self.topology, self.coords)

    def remove_solvent_and_ions(self):
        """ Remove the solvent and ions from this protein.
        """
        atoms = []
        for at in self.topology.iter_atoms():
            if at.residue.name not in SOLVENT_AND_IONS:
                atoms.append(at.index)

        self._topology = self._topology.subset(atoms)
        self._coords = self.coords[:, atoms, :]
        self._validate_coords(self.topology, self.coords)

    @staticmethod
    def _validate_coords(topology, coords):
        """ Verifies that the coordinates array has correct shape
            and that number of atoms match with the topology
        """
        if len(coords.shape) != 3:
            raise ValueError(f"Incorrect shape {coords.shape}")

        if coords.shape[1] != topology.n_atoms:
            raise ValueError(f"Incorrect number of atoms {coords.shape[1]}. "
                             f"Topology has {topology.n_atoms} atoms")
