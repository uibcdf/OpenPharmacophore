import abc
import numpy as np
import pyunitwizard as puw
from typing import List
from openpharmacophore.molecular_systems.chem_feats import ChemFeatContainer, mol_chem_feats
from openpharmacophore.molecular_systems.convert import topology_to_mol, ligand_to_topology
from openpharmacophore.molecular_systems.hbonds import protein_ligand_hbonds
from openpharmacophore.utils import maths


class AbstractBindingSite(abc.ABC):

    @abc.abstractmethod
    def get_atoms(self, frame: int) -> np.ndarray:
        raise NotImplementedError

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
        protein : openpharmacophore.Protein

        ligand : openpharmacophore.Ligand
    """

    def __init__(self, protein, ligand):
        self._protein = protein
        self._ligand = ligand

        self._bsite = None  # type: openpharmacophore.Protein
        self._bsite_mol = None  # type: rdkit.Chem.Mol
        self._has_ligand = None

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
        aminoacids, lig_res_ind = self.get_residues(frame)
        if self._ligand_has_hyd(lig_res_ind):
            # In this case we do not need to concatenate the ligand
            aminoacids.append(lig_res_ind)
            atoms = self._protein.topology.get_residues_atoms(aminoacids)
            self._has_ligand = True
        elif not self._ligand_has_hyd(lig_res_ind) or lig_res_ind is None:
            # We extract the binding site topology and concatenate the ligand topology
            atoms = self._protein.topology.get_residues_atoms(aminoacids)
            self._has_ligand = False
        else:
            assert False, "This should not happen."

        self._bsite = self._protein.slice(atoms, frame)

    def get_chem_feats(self, frame):
        """ Get the chemical features of the binding site at the specified frame.

            Parameters
            ----------
            frame : int

            Returns
            -------
            ChemFeatContainer
        """
        self._get_binding_site(frame)
        self._bsite_mol = topology_to_mol(
            self._bsite.topology,
            self._bsite.coords,
            remove_hyd=True)
        return mol_chem_feats(
            self._bsite_mol, self._bsite.coords[0], feat_def="protein"
        )

    def get_hydrogen_bonds(self, frame):
        """ Get the hydrogen bonds between the protein and the ligand.

            Parameters
            ----------
            frame : int

            Returns
            --------
            ChemFeatContainer
        """
        have_hyd = self._ligand.has_hydrogens and self._protein.has_hydrogens
        h_bonds = ChemFeatContainer()
        if have_hyd:
            if not self._has_ligand:
                # Add the ligand to the bsite, so we can get hydrogen bonds
                ligand_top = ligand_to_topology(self._ligand.to_rdkit())
                ligand_coords = np.expand_dims(self._ligand.get_conformer(frame), axis=0)
                self._bsite.concatenate(ligand_top, ligand_coords)

            h_bonds.add_feats(protein_ligand_hbonds(self._bsite))
        return h_bonds

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


class BindingSite(AbstractBindingSite):
    """ Class to extract and find chemical features in the binding
        site of a protein.
    """
    pass
