from rdkit.Chem import AllChem as Chem
import pyunitwizard as puw
import mdtraj as mdt
import tempfile

from openpharmacophore.molecular_systems.exceptions import DifferentNumAtomsError


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

    def fix_bond_order(self, smiles):
        """ Fix the bond order of a ligand. This is necessary when a ligand
            is extracted from a Protein object.

            Parameters
            ----------
            smiles : str
                Smiles of the ligand
        """
        template = Chem.MolFromSmiles(smiles)
        if self._mol.GetNumAtoms() != template.GetNumAtoms():
            # Removing Hs involved in double bonds can help with the
            # matching
            template = Chem.RemoveAllHs(template)
            if self.n_atoms != template.GetNumAtoms():
                raise DifferentNumAtomsError(
                    self.n_atoms, template.GetNumAtoms()
                )

        self._mol = Chem.AssignBondOrdersFromTemplate(template, self._mol)

    def has_aromatic_bonds(self):
        """ Returns true if the ligand has aromatic bonds.

            Returns
            -------
            bool
        """
        return any(
            [b.GetBondTypeAsDouble() == 1.5 for b in self._mol.GetBonds()]
        )

    def has_double_bonds(self):
        """ Returns true if the ligand has double bonds.

           Returns
           -------
           bool
       """
        return any(
            [b.GetBondTypeAsDouble() == 2.0 for b in self._mol.GetBonds()]
        )


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
        xyz=puw.get_value(coords, "nanometers")[0],
        topology=topology.top
    )
    # TODO: create ligand without using files
    pdb_file = tempfile.NamedTemporaryFile()
    traj.save_pdb(pdb_file.name)
    pdb_file.seek(0)

    mol = Chem.MolFromPDBFile(pdb_file.name)
    pdb_file.close()
    assert mol is not None, "Failed to create ligand"

    ligand = Ligand(mol)
    if coords.shape[0] > 1:
        ligand.add_conformers(coords)
    return ligand


def create_ligand_set(filename: str):
    raise NotImplementedError


def smiles_from_pdb_id(pdb_id: str):
    raise NotImplementedError
