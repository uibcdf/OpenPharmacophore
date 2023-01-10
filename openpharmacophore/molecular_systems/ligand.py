import mdtraj as mdt
import pyunitwizard as puw
from rdkit.Chem import AllChem as Chem
from pathlib import Path
import pickle
import tempfile

from openpharmacophore.molecular_systems.exceptions import DifferentNumAtomsError, LigandCreationError


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

    @classmethod
    def from_string(cls, string, form):
        """ Create a ligand from a smiles, smarts, inchi, pdb block,
            or mol2 block


            Parameters
            ----------
            string: str
                The string that encodes the ligand

            form : str
                Can be "smi", "smarts", "inchi", "mol2", "pdb".

            Returns
            -------
            Ligand

        """
        if form == "smi":
            mol = Chem.MolFromSmiles(string)
        elif form == "smarts":
            mol = Chem.MolFromSmarts(string)
        elif form == "inchi":
            mol = Chem.MolFromInchi(string)
        elif form == "mol2":
            mol = Chem.MolFromMol2Block(string)
        elif form == "pdb":
            mol = Chem.MolFromPDBBlock(string)
        else:
            raise NotImplementedError(f"Cannot create a ligand from form {form}")

        if mol is None:
            raise LigandCreationError(string, form)
        return cls(mol)

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

    def has_hydrogens(self):
        """ Returns true if the ligand has any hydrogen atoms

            Returns
            -------
            bool
        """
        return any([
            a.GetSymbol() == "H" for a in self._mol.GetAtoms()
        ])

    def add_hydrogens(self):
        """ Add hydrogens to this ligand. Modifies the ligand in place.
        """
        self._mol = Chem.AddHs(self._mol, addCoords=True, addResidueInfo=True)


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


def _load_pdb_id_mapper():
    """ Returns a dictionary that maps pdb ids to smiles.

        Returns
        -------
        dict : [str, str]

    """
    path = Path(__file__).parent / "pdb_to_smi.pickle"
    with open(path, "rb") as fp:
        mapper = pickle.load(fp)
    return mapper


def smiles_from_pdb_id(pdb_id: str, mapper=None):
    """ Obtain a smiles from a ligand pdb id.

        Parameters
        ----------
        pdb_id : str
            The ligand id with oir without the chain i.e. the following are
            equivalent "DAO" and "DAO:A"

        mapper : dict[str, str], optional
            Maps the pdb id to the smiles. If None is provided a default
            mapper is loaded by openpharmacophore.

        Returns
        -------
        str

    """
    if mapper is None:
        mapper = _load_pdb_id_mapper()
    pdb_id = pdb_id.split(":")[0]
    return mapper[pdb_id]
