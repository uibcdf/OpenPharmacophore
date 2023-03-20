import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Geometry.rdGeometry import Point3D
import pyunitwizard as puw

from copy import deepcopy
from pathlib import Path
import pickle

from openpharmacophore.molecular_systems.exceptions import DifferentNumAtomsError, LigandCreationError
from openpharmacophore.molecular_systems.convert import topology_to_mol
from openpharmacophore.molecular_systems.chem_feats import get_indices, mol_chem_feats, \
    aromatic_chem_feats, donor_chem_feats
from openpharmacophore.molecular_systems.chem_feats import ChemFeatContainer, SMARTS_LIGAND
from openpharmacophore.utils.conformers import ConformerGenerator


class LigandWithHsError(ValueError):
    """ Exception raised when trying to add Hs to a ligand that
        already contains hydrogen.
    """
    pass


class Ligand:
    """ Represents a ligand. Wrapper for rdkit Mol object
    """
    def __init__(self, mol):
        self._mol = mol  # type: Chem.Mol

        self._feat_ind = {}  # type: dict[str, list[tuple]]
        self._lig_id = None
        self._has_hyd = Ligand._mol_has_hyd(self._mol)

    @property
    def n_atoms(self) -> int:
        return self._mol.GetNumAtoms()

    @property
    def n_heavy_atoms(self) -> int:
        return self._mol.GetNumHeavyAtoms()

    @property
    def n_conformers(self) -> int:
        return self._mol.GetNumConformers()

    @property
    def n_bonds(self) -> int:
        return self._mol.GetNumBonds()

    @property
    def has_hydrogens(self) -> bool:
        """ Returns true if the ligand has any hydrogen atoms

            Returns
            -------
            bool
        """
        return self._has_hyd

    @property
    def has_conformers(self) -> bool:
        return self._mol.GetNumConformers() > 0

    @property
    def lig_id(self) -> str:
        """ ID of the ligand if it was extracted from a protein.
        """
        return self._lig_id

    @lig_id.setter
    def lig_id(self, ligand_id: str):
        self._lig_id = ligand_id

    @property
    def feat_ind(self):
        return self._feat_ind

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

        coords = puw.get_value(coords, "angstroms")
        n_conformers = coords.shape[0]
        n_atoms = coords.shape[1]
        for ii in range(n_conformers):
            conformer = Chem.Conformer(n_atoms)
            for jj in range(n_atoms):
                pos = [float(co) for co in coords[ii, jj, :]]
                conformer.SetAtomPosition(jj, Point3D(*pos))
            self._mol.AddConformer(conformer, assignId=True)

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
        return puw.quantity(
            self._mol.GetConformer(conf_ind).GetPositions(),
            "angstroms",
        )

    def to_rdkit(self):
        """ Returns the ligand as an rdkit molecule.

            Returns
            -------
            rdkit.Chem.Mol
        """
        return self._mol

    def get_chem_feats(self, conf_ind, types=None, indices=False):
        """ Find chemical features in this ligand.

            Parameters
            ----------
            conf_ind : int
                Index of the conformer.

            types : set[str], optional
                The chemical features that will be searched for.

            indices : bool, default=False
                Whether to store the indices of the atoms in the chemical features.

            Returns
            -------
            ChemFeatContainer
        """
        if len(self._feat_ind) == 0:
            self._feat_ind = get_indices(self._mol, feat_def=SMARTS_LIGAND, types=types)
        return mol_chem_feats(self._feat_ind, self.get_conformer(conf_ind), store_indices=indices)

    def get_chem_feats_with_directionality(self, conf_ind, types=None):
        """ Find chemical features in this ligand. Hydrogen bond donors obtained
            with this method will include information about hydrogen atoms, and
            aromatic ring will contain their normal vector.

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
        if len(self._feat_ind) == 0:
            self._feat_ind = get_indices(self._mol, feat_def=SMARTS_LIGAND, types=types)

        # Exclude donor and aromatics and process them separately
        donors = self._feat_ind.pop("hb donor")
        aromatics = self._feat_ind.pop("aromatic ring")

        conformer = self.get_conformer(conf_ind)
        container = mol_chem_feats(self._feat_ind, conformer)
        container.add_feats(donor_chem_feats(donors, conformer, self._mol))
        container.add_feats(aromatic_chem_feats(aromatics, conformer))

        self._feat_ind["hb donor"] = donors
        self._feat_ind["aromatic ring"] = aromatics

        return container

    def draw(self):
        """ Draw the ligand.
        """
        mol_copy = deepcopy(self._mol)
        # Remove conformers for pretty drawing
        mol_copy.RemoveAllConformers()
        return mol_copy

    def generate_conformers(self, n_confs):
        """ Generate conformers for this ligand. Removes current
            conformers if any
        """
        self._mol.RemoveAllConformers()
        conf_gen = ConformerGenerator(n_confs)
        self._mol = conf_gen(self._mol)

    @staticmethod
    def _mol_has_hyd(mol):
        """ Returns true if the molecule contains hydrogens.
        """
        if any(atom.GetSymbol() == "H" for atom in mol.GetAtoms()):
            n_implicit_hs = sum(a.GetNumImplicitHs() for a in mol.GetAtoms())
            n_actual_hs = sum(a.GetSymbol() == "H" for a in mol.GetAtoms())
            return not (n_actual_hs < n_implicit_hs)
        return False


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
    assert len(coords.shape) == 3, f"Coordinates array has incorrect shape {coords.shape}"
    mol = topology_to_mol(
        topology, puw.get_value(coords, "nanometers")[0], remove_hyd)
    ligand = Ligand(mol)
    ligand.lig_id = topology.get_residue(0).name

    if coords.shape[0] > 1:
        ligand.add_conformers(coords[1:, :, :])

    return ligand


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
