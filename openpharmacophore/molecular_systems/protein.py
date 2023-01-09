from openpharmacophore.molecular_systems import Topology


class Protein:
    """ Represents a protein and optionally its different structures
        if it comes from a MD trajectory.

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
