from rdkit.Chem import AllChem as Chem
import pyunitwizard as puw
import mdtraj as mdt
import tempfile


class Ligand:
    """ Represents a ligand. Wrapper for rdkit Mol object
    """
    def __init__(self, mol):
        self._mol = mol  # type: Chem.Mol

    @property
    def n_atoms(self) -> int:
        return self._mol.GetNumAtoms()

    @property
    def n_conformers(self) -> int:
        return self._mol.GetNumConformers()

    @property
    def n_bonds(self) -> int:
        return self._mol.GetNumBonds()


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
    traj = mdt.Trajectory(
        xyz=puw.get_value(coords, "nanometers"),
        topology=topology.top
    )
    # TODO: create ligand without using files
    pdb_file = tempfile.NamedTemporaryFile()
    traj.save_pdb(pdb_file.name)
    pdb_file.seek(0)

    mol = Chem.MolFromPDBFile(pdb_file.name)
    assert mol is not None, "Failed to create ligand"

    pdb_file.close()

    return Ligand(mol)


def create_ligand_set(filename: str):
    pass
