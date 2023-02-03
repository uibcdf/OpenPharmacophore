import abc
import numpy as np
import pyunitwizard as puw
from typing import List
import openpharmacophore.molecular_systems.chem_feats as cf
from openpharmacophore.molecular_systems.convert import topology_to_mol
from openpharmacophore.utils import maths


class AbstractBindingSite(abc.ABC):

    @abc.abstractmethod
    def get_atoms(self, frame: int) -> np.ndarray:
        raise NotImplementedError

    @abc.abstractmethod
    def get_residues(self, frame: int) -> List[int]:
        raise NotImplementedError

    @abc.abstractmethod
    def get_chem_feats(self, frame: int) -> cf.ChemFeatContainer:
        raise NotImplementedError

    @abc.abstractmethod
    def to_rdkit(self):
        raise NotImplementedError


BS_DIST_MAX = puw.quantity(0.85, "nanometers")


class ComplexBindingSite(AbstractBindingSite):
    """ Class to extract and find chemical features in the binding
        site of a protein-ligand complex.

        Parameters
        ----------
        protein : openpharmacophore.Protein

        ligand : openpharmacophore.Ligand
    """

    def __init__(self, protein, ligand):
        self._protein = protein
        self._ligand = ligand

        self._bsite_mol = None  # type: rdkit.Chem.Mol

    def get_atoms(self, frame):
        """ Get the indices of the atoms in the binding site. This will include the ligand
            atoms as well.

            Parameters
            ----------
            frame : int

            Returns
            -------
            np.ndarray
                Array with the atom indices
        """
        conformer = self._ligand.get_conformer(frame)
        lig_centroid = self._ligand_centroid(conformer)
        max_extent = self._ligand_max_extent(conformer, lig_centroid)
        bs_cutoff = max_extent + BS_DIST_MAX

        return self._protein.atoms_at_distance(
            frame, lig_centroid, max_dist=bs_cutoff, min_dist=0)

    def get_residues(self, frame):
        """ Get the indices of the residues in the protein that correspond
            to the binding site.

            The binding site is defined as the sphere with centroid at the ligand
            centroid with a radius of ligand maximum extent plus a cutoff distance.

            Parameters
            ----------
            frame : int

            Returns
            -------
            list[int]
                Indices of the residues that encompass the binding site
        """
        atoms = self.get_atoms(frame)
        return self._protein.topology.get_bs_residues(atoms)

    def _get_binding_site(self, frame):
        """ Create a Protein object representing the binding site

            Parameters
            ----------
            frame : int

            Returns
            -------
            openPharmacophore.Protein
                The binding site
        """
        aminoacids, _ = self.get_residues(frame)
        atoms = self._protein.topology.get_residues_atoms(aminoacids)
        return self._protein.slice(atoms, frame)

    def get_chem_feats(self, frame, types=None):
        """ Get the chemical features of the binding site at the specified frame.

            Parameters
            ----------
            frame : int

            types : set[str], optional
                Types of features to search for. If none all feat types
                will be used.

            Returns
            -------
            ChemFeatContainer
        """
        b_site = self._get_binding_site(frame)
        self._bsite_mol = topology_to_mol(
            b_site.topology, b_site.coords, remove_hyd=True)

        indices = cf.get_indices(self._bsite_mol, feat_def=cf.SMARTS_PROTEIN, types=types)
        # Donors and aromatics need to be processed differently because they require extra data
        donors = indices.pop("hb donor")
        aromatics = indices.pop("aromatic ring")

        container = cf.mol_chem_feats(indices, b_site.coords[0])
        container.add_feats(cf.donor_chem_feats(donors, b_site.coords[0], self._bsite_mol))
        container.add_feats(cf.aromatic_chem_feats(aromatics, b_site.coords[0]))

        return container

    @staticmethod
    def _ligand_centroid(coords):
        """ Get the centroid of the ligand.

            Parameters
            ---------
            coords : puw.Quantity
                Shape (n_atoms, 3)
            Returns
            -------
            puw.Quantity
                Shape (3, 3)
        """
        return np.mean(coords, axis=0)

    @staticmethod
    def _ligand_max_extent(coords, centroid):
        """ Computes the maximum distance from the ligand centroid
            to any of its atoms.

            Parameters
            ----------
            coords : puw.Quantity
                Shape (n_atoms, 3)

            centroid : puw.Quantity
                Shape (3, 3)

            Returns
            -------
            puw.Quantity
                A scalar

        """
        return maths.maximum_distance(centroid, coords)

    def _ligand_has_hyd(self, lig_res_ind):
        """ Returns true of the ligand has hydrogens.
        """
        return lig_res_ind is not None and \
            self._protein.topology.residue_has_hyd(lig_res_ind)

    def to_rdkit(self):
        """ Returns the binding site as an rdkit molecule.

            Returns
            -------
            rdkit.Chem.Mol
        """
        return self._bsite_mol


class BindingSite(AbstractBindingSite):
    """ Class to extract and find chemical features in the binding
        site of a protein.
    """
    pass
