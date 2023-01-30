from openpharmacophore import PharmacophoricPoint, Pharmacophore
from openpharmacophore.utils import maths
import networkx as nx
import numpy as np
import pyunitwizard as puw
import json
import re
import requests
import tempfile


class LigandReceptorPharmacophore:
    """ Class to store, and extract pharmacophores from protein-ligand complexes.

        The pharmacophores can be extracted from a PDB structure or from a molecular
        dynamics trajectory.

    """
    # Values from ligandscout and plip
    HYD_DIST_MAX = puw.quantity(0.5, "nanometers")
    HYD_MERGE_DIST = puw.quantity(0.2, "nanometers")  # value from pharmer

    CHARGE_DIST_MAX = puw.quantity(0.56, "nanometers")

    PISTACK_DIST_MAX = puw.quantity(0.75, "nanometers")
    PISTACK_OFFSET_MAX = puw.quantity(0.20, "nanometers")
    PISTACK_ANG_DEV = 30  # degrees

    def __init__(self):
        self._pharmacophores = []  # type: list[Pharmacophore]
        self._pl_complex = None

    @property
    def num_frames(self):
        return len(self._pharmacophores)

    def extract(self):
        raise NotImplementedError

    def _hbond_donor_pharmacophoric_points(self, h_bonds, frame):
        """ Compute hydrogen bond donor pharmacophoric points from
            protein-ligand interactions.

            Parameters
            -----------
            h_bonds : np.ndarray
                Array with the indices of the atoms involved in hydrogen bond.
                Each row contains three integer indices, (d_i, h_i, a_i), such
                that d_i is the index of the donor atom, h_i the index of the
                hydrogen atom, and a_i the index of the acceptor atom.
                Shape = (n_h_bonds, 3)

            frame : int
                The frame of the trajectory.
        """
        radius = puw.quantity(1.0, "angstroms")
        for bond in h_bonds:
            # There is an acceptor in the receptor
            if bond[0] in self._pl_complex.lig_indices:
                direction = puw.get_value(self._pl_complex.coords[frame, bond[2], :] -
                                          self._pl_complex.coords[frame, bond[0], :])
                pharma_point = PharmacophoricPoint(
                    "hb donor", self._pl_complex.coords[frame, bond[0], :],
                    radius, direction
                )
                self._pharmacophores[frame].add(pharma_point)

    def _hbond_acceptor_pharmacophoric_points(self, h_bonds, frame):
        """ Compute hydrogen bond acceptor pharmacophoric points from
            protein-ligand interactions.

            Parameters
            -----------
            h_bonds : np.ndarray
                Array with the indices of the atoms involved in hydrogen bond.
                Each row contains three integer indices, (d_i, h_i, a_i), such
                that d_i is the index of the donor atom, h_i the index of the
                hydrogen atom, and a_i the index of the acceptor atom.
                Shape = (n_h_bonds, 3)

            frame : int
                The frame of the trajectory.
        """
        radius = puw.quantity(1.0, "angstroms")
        for bond in h_bonds:
            if bond[2] in self._pl_complex.lig_indices:
                # There is a donor in the ligand
                direction = puw.get_value(self._pl_complex.coords[frame, bond[0], :] -
                                          self._pl_complex.coords[frame, bond[2], :])
                pharma_point = PharmacophoricPoint(
                    "hb acceptor", self._pl_complex.coords[frame, bond[2], :],
                    radius, direction
                )
                self._pharmacophores[frame].add(pharma_point)

    def _aromatic_pharmacophoric_points(self, lig_centers, lig_indices,
                                        rec_centers, rec_indices, frame):
        """ Compute aromatic pharmacophoric points from
            protein-ligand interactions.

            Parameters
            -----------
            lig_centers : list[puw.Quantity]
                Centroids of the aromatic rings in the ligand.

            lig_indices : list[list[int]]
                Indices of the rings atoms in the ligand.

            rec_centers : list[puw.Quantity]
                Centroids of aromatic rings in the receptor.

            rec_indices : list[list[int]]
                Indices of the rings atoms in the receptor.

            frame : int
                The frame of the trajectory.

        """
        # Calculate pistack interactions and create pharmacophoric points
        radius = puw.quantity(1.0, "angstroms")
        for ii in range(len(lig_centers)):
            for jj in range(len(rec_centers)):
                if maths.points_distance(rec_centers[jj], lig_centers[ii]) <= self.PISTACK_DIST_MAX:
                    # Calculate deviation from ideal angle by taking the angle between the normals
                    # defined by the planes of each ring
                    lig_normal = maths.ring_normal(lig_indices[ii], self._pl_complex.coords[frame], lig_centers[ii])
                    rec_normal = maths.ring_normal(rec_indices[jj], self._pl_complex.coords[frame], rec_centers[jj])
                    angle = maths.angle_between_normals(lig_normal, rec_normal)
                    assert 0 <= angle <= 360, f"Angle is {angle}"

                    if 0 <= angle <= self.PISTACK_ANG_DEV or \
                            90 - self.PISTACK_ANG_DEV <= angle <= 90 + self.PISTACK_ANG_DEV:

                        # Project ring centers into the other plane and calculate offset
                        rec_proj = maths.point_projection(lig_normal, lig_centers[ii], rec_centers[jj])
                        lig_proj = maths.point_projection(rec_normal, rec_centers[jj], lig_centers[ii])
                        offset = min(maths.points_distance(lig_proj, rec_centers[jj]),
                                     maths.points_distance(rec_proj, lig_centers[ii]))
                        if offset <= self.PISTACK_OFFSET_MAX:
                            direction = puw.get_value(rec_centers[jj] - lig_centers[ii])
                            pharma_point = PharmacophoricPoint(
                                "aromatic ring", lig_centers[ii], radius, direction)
                            self._pharmacophores[frame].add(pharma_point)

    def _hydrophobic_pharmacophoric_points(self, ligand_centers, receptor_centers, frame):
        """ Compute hydrophobic pharmacophoric points from protein-ligand interactions.

            Typically, there will be a lot of hydrophobic points so, they are passed
            to a merging procedure.

            Parameters
            -----------
            ligand_centers : list[puw.Quantity]
                Centroids of the hydrophobic areas in the ligand.

            receptor_centers : list[puw.Quantity]
                Centroids of hydrophobic areas in the receptor.

            frame : int
                The frame of the trajectory.

        """
        centers = []
        radius = puw.quantity(1.0, "angstroms")
        for lig_center in ligand_centers:
            for prot_center in receptor_centers:
                dist = maths.points_distance(lig_center, prot_center)
                if dist < self.HYD_DIST_MAX:
                    centers.append(lig_center)

        points_clustered = self._merge_hydrophobic_points(centers, radius)
        for p in points_clustered:
            self._pharmacophores[frame].add(p)

    @staticmethod
    def _merge_hydrophobic_points(centers, radius):
        """ Merge group of hydrophobic points close to each other.

            Parameters
            ----------
            centers : list[puw.Quantity]
            radius : puw.Quantity

            Returns
            -------
            list[PharmacophoricPoint]
        """
        # Create a graph of the hydrophobic features, were each node represents
        # a feature and an edge is added between two nodes only if their distance
        # is < HYD_MERGE_DIST.
        hyd_graph = nx.Graph()
        for ii in range(len(centers)):
            hyd_graph.add_node(ii)

        for ii in range(len(centers)):
            for jj in range(ii + 1, len(centers)):
                dist = maths.points_distance(centers[ii], centers[jj])
                if dist <= LigandReceptorPharmacophore.HYD_MERGE_DIST:
                    hyd_graph.add_edge(ii, jj)

        # Find each maximum clique and group all nodes within a clique
        cliques_iter = nx.find_cliques(hyd_graph)
        points = []
        for clique in cliques_iter:
            clique_centers = [centers[ii] for ii in clique]
            center = np.mean(np.stack(clique_centers), axis=0)
            pharma_point = PharmacophoricPoint("hydrophobicity", center, radius)
            points.append(pharma_point)
        return points

    def _charge_pharmacophoric_points(self, ligand_centers, receptor_centers,
                                      charge_type, frame):
        """ Compute positive or negative charge pharmacophoric points from
            protein-ligand interactions.

            Parameters
            -----------
            ligand_centers : list[puw.Quantity]
                Centroids of the positive or negative areas in the ligand.

            receptor_centers : list[puw.Quantity]
                Centroids of areas with opposite sign in the receptor.

            frame : int
                The frame of the trajectory.

        """
        radius = puw.quantity(1.0, "angstroms")
        for lig_center in ligand_centers:
            for prot_center in receptor_centers:
                if maths.points_distance(lig_center, prot_center) < self.CHARGE_DIST_MAX:
                    pharma_point = PharmacophoricPoint(charge_type, lig_center, radius)
                    self._pharmacophores[frame].add(pharma_point)

    def __len__(self):
        return len(self._pharmacophores)

    def __getitem__(self, frame):
        """
            Returns
            -------
            Pharmacophore
        """
        return self._pharmacophores[frame]
