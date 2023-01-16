import mdtraj as mdt
from rdkit import Chem
import tempfile


def topology_to_mol(topology, coords, remove_hyd=True):
    """ Create a ligand from a topology object and an array of coordinates

           Parameters
           ----------
           topology : Topology
               Topology of the ligand.

           coords: np.ndarray
               Array of shape (n_atoms, 3).

           remove_hyd : bool
               Whether to remove the hydrogens from the molecule.

           Returns
           -------
           mol: rdkit.Chem.Mol

   """
    traj = mdt.Trajectory(
        xyz=coords,
        topology=topology.top
    )
    # TODO: create molecule without using files
    pdb_file = tempfile.NamedTemporaryFile()
    traj.save_pdb(pdb_file.name)
    pdb_file.seek(0)

    mol = Chem.MolFromPDBFile(pdb_file.name, removeHs=remove_hyd)
    pdb_file.close()

    assert mol is not None, "Failed to create molecule"
    mol.RemoveAllConformers()
    return mol
