
class Ligand:
    pass


class LigandSet:
    pass


def ligand_from_topology(topology, coords):
    """ Create a ligand from a topology object and an array of coordinates

        Parameters
        ----------
        topology : Topology
            Topology of the ligand.

        coords: Quantity
            Array of shape (n_conformers, n_atoms, 3).

        Returns
        -------
        Ligand

    """
    pass


def create_ligand_set(filename: str):
    pass
