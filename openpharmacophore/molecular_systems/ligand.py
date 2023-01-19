import numpy as np
from rdkit.Chem import AllChem as Chem
import pyunitwizard as puw
from pathlib import Path
import pickle

from openpharmacophore.molecular_systems.exceptions import DifferentNumAtomsError, LigandCreationError
from openpharmacophore.molecular_systems.convert import topology_to_mol
from openpharmacophore.molecular_systems.chem_feats import get_indices, mol_chem_feats
from openpharmacophore.molecular_systems.chem_feats import ChemFeatContainer, SMARTS_LIGAND


class LigandWithHsError(ValueError):
    """ Exception raised when trying to add Hs to a ligand that
        already contains hydrogen.
    """
    pass


class Ligand:
    """ Represents a ligand. Wrapper for rdkit Mol object
    """
    def __init__(self, mol, conformer_coords=None):
        self._mol = mol  # type: Chem.Mol
        # quantity of shape (n_conformers, n_atoms, 3)
        self._conformers = conformer_coords

        self._feat_ind = None
        self._has_hyd = None
        self._lig_id = None

    @property
    def n_atoms(self) -> int:
        return self._mol.GetNumAtoms()

    @property
    def n_heavy_atoms(self) -> int:
        return self._mol.GetNumHeavyAtoms()

    @property
    def n_conformers(self) -> int:
        if not self.has_conformers:
            return 0
        return self._conformers.shape[0]

    @property
    def n_bonds(self) -> int:
        return self._mol.GetNumBonds()

    @property
    def has_hydrogens(self):
        """ Returns true if the ligand has any hydrogen atoms

            Returns
            -------
            bool
        """
        if self._has_hyd is None:
            self._has_hyd = any([
                a.GetSymbol() == "H" for a in self._mol.GetAtoms()
            ])
        return self._has_hyd

    @property
    def has_conformers(self):
        return self._conformers is not None

    @property
    def lig_id(self) -> str:
        """ ID of the ligand if it was extracted from a protein.
        """
        return self._lig_id

    @lig_id.setter
    def lig_id(self, ligand_id: str):
        self._lig_id = ligand_id

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
        if not self._has_hyd and self.n_atoms != template.GetNumAtoms():
            # Removing Hs involved in double bonds can help with the
            # matching
            template = Chem.RemoveAllHs(template)
            if self.n_atoms != template.GetNumAtoms():
                raise DifferentNumAtomsError(
                    self.n_atoms, template.GetNumAtoms()
                )
        elif self._has_hyd:
            # TODO: can we do the matching without adding Hydrogens to speed the process?
            template = Chem.AddHs(template)

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

    def add_hydrogens(self):
        """ Add hydrogens to this ligand. Modifies the ligand in place.

            Can only add hydrogens to a ligand that only has one conformer.
        """
        if self.has_hydrogens:
            raise LigandWithHsError("This ligand already contains hydrogens")

        self._mol = Chem.AddHs(self._mol, addCoords=True, addResidueInfo=True)
        self._has_hyd = True

        if self.has_conformers and self.n_conformers == self._mol.GetNumConformers():
            self._conformers = puw.quantity(
                self._get_molecule_coordinates(self._mol), "angstroms"
            )

        if self.n_conformers != self._mol.GetNumConformers():
            # When the rdkit mol has no conformers and the conformers array is not null
            # the hydrogens will not have coordinates.
            raise ValueError("Failed to add Hs")

        # We do not need to store the conformers in the rdkit molecule anymore
        self._mol.RemoveAllConformers()

    @staticmethod
    def _get_molecule_coordinates(molecule):
        """ Returns an array with the positions of the atoms in each
            conformer of the molecule.

            Parameters
            ----------
            molecule : rdkit.Chem.Mol

            Returns
            -------
            np.ndarray
                Array of shape (n_conformers, n_atoms, 3)
        """
        conformer: Chem.Conformer

        n_conformers = molecule.GetNumConformers()
        n_atoms = molecule.GetNumAtoms()
        new_coords = np.zeros((n_conformers, n_atoms, 3))
        for ii in range(n_conformers):
            conformer = molecule.GetConformer(ii)
            new_coords[ii] = conformer.GetPositions()

        return new_coords

    def add_conformers(self, coords):
        """ Add conformers to a ligand.

            Parameters
            ----------
            coords : puw.Quantity
                A quantity of shape (n_conformers, n_atoms, 3)

        """
        if len(coords.shape) != 3 or coords.shape[1] != self.n_atoms:
            shape = ("x", self.n_atoms, 3)
            raise ValueError(f"Incorrect shape {coords.shape}"
                             f"Expected shape {shape}")

        if self._conformers is None:
            self._conformers = coords
        else:
            self._conformers = np.concatenate((self._conformers, coords))

    def get_conformer(self, conf_ind):
        """ Get the coordinates of the specified conformer

            Parameters
            ----------
            conf_ind : int
                Index of the conformer

            Returns
            -------
            puw.Quantity
                A quantity of shape (n_atoms, 3)
        """
        return self._conformers[conf_ind]

    def to_rdkit(self):
        """ Returns the ligand as an rdkit molecule.

            Returns
            -------
            rdkit.Chem.Mol
        """
        return self._mol

    def get_chem_feats(self, conf_ind, types=None):
        """ Find chemical features in this ligand.

            Parameters
            ----------
            conf_ind : int
                Index of the conformer.

            types : set[str], optional
                The chemical features that will be searched for.

            Returns
            -------
            ChemFeatContainer
        """
        if self._feat_ind is None:
            self._feat_ind = get_indices(self._mol, feat_def=SMARTS_LIGAND, types=types)
        return mol_chem_feats(self._feat_ind, self.get_conformer(conf_ind))


class LigandSet:
    pass


def ligand_from_topology(topology, coords, remove_hyd=True):
    """ Create a ligand from a topology object and an array of coordinates

        Parameters
        ----------
        topology : Topology
            Topology of the ligand.

        coords: Quantity
            Array of shape (n_conformers, n_atoms, 3).

        remove_hyd : bool
            Whether to remove the hydrogens from the ligand.

        Returns
        -------
        Ligand

    """
    assert topology.n_residues == 1, f"Topology has {topology.n_residues}"
    mol = topology_to_mol(topology,
                          puw.get_value(coords, "nanometers")[0],
                          remove_hyd)
    if remove_hyd:
        non_hyd = topology.non_hyd_indices()
        coords = coords[:, non_hyd, :]

    ligand = Ligand(mol, coords)
    ligand.lig_id = topology.get_residue(0).name
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
