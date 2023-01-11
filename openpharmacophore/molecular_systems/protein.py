from openpharmacophore.molecular_systems.topology import Topology
from openpharmacophore.molecular_systems.ligand import ligand_from_topology
from openpharmacophore.utils.maths import delete
from openmm.app import Modeller
import pyunitwizard as puw
import numpy as np


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
        self._topology = topology
        self._coords = coords

    @property
    def n_atoms(self):
        return self._topology.n_atoms

    @property
    def has_hydrogens(self) -> bool:
        return self._topology.has_hydrogens()

    @property
    def has_ligands(self) -> bool:
        return self._topology.has_ligands()

    @property
    def ligand_ids(self):
        return self._topology.ligand_ids()

    def get_ligand(self, ligand_id, remove=False, remove_hyd=True):
        """ Extract a ligand assuming there is one.

            Parameters
            ----------
            ligand_id : str
                ID of the ligand. ID is composed by pdb id followed by chain name,
                i.e. "EST:B".

            remove : bool, default=False
                Whether to remove the ligand from the topology.

            remove_hyd : bool
                Whether to remove teh hydrogens from the ligand.

            Returns
            -------
            Ligand
                A ligand. Its bond orders may be incorrect.
        """
        lig_name, chain = ligand_id.split(":")
        lig_indices = self._topology.get_residue_indices(lig_name, chain)
        assert len(lig_indices) > 0, f"Ligand {ligand_id} could not be extracted"

        lig_coords = self._coords[:, lig_indices, :]
        ligand_top = self._topology.subset(lig_indices)
        ligand = ligand_from_topology(ligand_top, lig_coords, remove_hyd)

        if remove:
            self._remove_ligand_by_indices(lig_indices)
        return ligand

    def _remove_ligand_by_indices(self, ligand_indices):
        """ Remove a ligand from this protein.

            Parameters
            ----------
            ligand_indices : list[int]

        """
        self._topology.remove_atoms(ligand_indices, inplace=True)
        self._coords = delete(self._coords, ligand_indices, axis=1)

    def _remove_ligand_by_id(self, ligand_id):
        """ Remove a ligand from this protein.

            Parameters
            ----------
            ligand_id : str
                ID of the ligand

        """
        lig_indices = self._topology.get_residue_indices(ligand_id)
        self._topology.remove_atoms(lig_indices, inplace=True)
        self._coords = delete(self._coords, lig_indices, axis=1)

    def add_hydrogens(self):
        """ Add hydrogens to the protein.
        """
        modeller = Modeller(
            topology=self._topology.to_openmm(),
            positions=puw.convert(self._coords[0],
                                  to_unit="nanometers",
                                  to_form="openmm.unit")
        )
        modeller.addHydrogens()
        self._topology, self._coords = modeller_to_topology(modeller)


def modeller_to_topology(modeller):
    """ Convert an openmm.Modeller to a mdtraj.Trajectory.

        Parameters
        ----------
        modeller : openmm.Modeller

        Returns
        -------
        Topology
            The topology

        coords : puw.Quantity
            Coordinates of the protein.
    """
    positions = modeller.getPositions()
    coords = puw.convert(positions, to_unit="angstroms", to_form="pint")
    coords = np.expand_dims(coords, axis=0)
    assert coords.shape == (1, len(positions), 3), f"Incorrect shape {coords.shape}"

    topology = Topology.from_openmm(modeller.getTopology())
    return topology, coords
