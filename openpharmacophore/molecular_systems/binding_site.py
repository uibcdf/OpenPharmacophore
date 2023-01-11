import abc
from typing import List
import numpy as np
import pyunitwizard as puw
from openpharmacophore.molecular_systems.chem_feats import ChemFeatContainer, residues_chem_feats
from openpharmacophore.utils import maths


class AbstractBindingSite(abc.ABC):

    @abc.abstractmethod
    def get_residues(self, frame: int) -> List[int]:
        raise NotImplementedError

    @abc.abstractmethod
    def get_chem_feats(self, frame: int) -> ChemFeatContainer:
        raise NotImplementedError


BS_DIST_MAX = puw.quantity(0.85, "nanometers")


class ComplexBindingSite(AbstractBindingSite):
    """ Class to extract and find chemical features in the binding
        site of a protein-ligand complex.

        Parameters
        ----------
        protein : Protein

        ligand : Ligand
    """
    def __init__(self, protein, ligand):
        self._protein = protein
        self._ligand = ligand

    def get_residues(self, frame):
        """ Get the indices of the residues in the protein that corresponds
            to the binding site.

            The binding site is defined as the sphere with centroid at the ligand
            centroid with a radius of ligand maximum extent plus a cutoff distance.


            Returns
            -------
            list[int]
                Indices of the residues that encompass the binding site
        """
        conformer = self._ligand.get_conformer(frame)
        lig_centroid = self._ligand_centroid(conformer)
        max_extent = self._ligand_max_extent(conformer, lig_centroid)
        bs_cutoff = max_extent + BS_DIST_MAX
        return self._protein.residues_at_distance(
            frame, lig_centroid, bs_cutoff, max_extent
        )

    def get_chem_feats(self, frame):
        """ Get the chemical features of the binding site at the specified frame.

            Parameters
            ----------
            frame : int

            Returns
            -------
            ChemFeatContainer
        """
        residues = self.get_residues(frame)
        return residues_chem_feats(residues)

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


class BindingSite(AbstractBindingSite):
    """ Class to extract and find chemical features in the binding
        site of a protein.
    """
    pass
