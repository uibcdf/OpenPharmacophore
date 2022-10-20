from .pharmacophore import Pharmacophore
from .pharmacophoric_point import PharmacophoricPoint
from .rdkit_pharmacophore import rdkit_pharmacophore
from .pl_complex import PLComplex
from ..io import (json_pharmacophoric_elements, ligandscout_xml_tree,
                  mol2_file_info, ph4_string)
from .._private_tools.exceptions import PDBFetchError
from ..utils import maths as maths
import networkx as nx
import numpy as np
import nglview as nv
import pyunitwizard as puw
import json
import re
import requests
import tempfile


class LigandReceptorPharmacophore(Pharmacophore):
    """ Class to store, and extract pharmacophores from protein-ligand complexes.

        The pharmacophores can be extracted from a pdb file or from a molecular
        dynamics simulation.

    """
    # Values from ligandscout and plip
    HB_DIST_MIN = puw.quantity(2.5, "angstroms")
    HB_DIST_MAX = puw.quantity(4.1, "angstroms")
    HB_DON_ANG_MIN = 100  # degrees

    HYD_DIST_MAX = puw.quantity(5.0, "angstroms")
    HYD_MERGE_DIST = puw.quantity(2.0, "angstroms")  # value from pharmer

    CHARGE_DIST_MAX = puw.quantity(5.6, "angstroms")

    PISTACK_DIST_MAX = puw.quantity(7.5, "angstroms")
    PISTACK_OFFSET_MAX = puw.quantity(2.0, "angstroms")
    PISTACK_ANG_DEV = 30  # degrees

    def __init__(self):
        # Pharmacophores will be stored as a list of pharmacophoric points.
        # A list for each pharmacophore
        self._pharmacophores = []
        self._pharmacophores_frames = []  # Contains the frame to which each pharmacophore belongs
        self._num_frames = 0

        self._pl_complex = None

    @property
    def num_frames(self):
        return self._num_frames

    @property
    def receptor(self):
        return self._pl_complex

    @staticmethod
    def _is_pdb_id(receptor):
        """ Check if the receptor is a PDB id.

            Parameters
            ----------
            receptor: str
                The receptor should be a string

            Returns
            -------
            bool
                Whether the receptor is a pdb id.
        """
        if len(receptor) == 4:
            pattern = re.compile('[0-9][a-zA-Z_0-9]{3}')
            if pattern.match(receptor):
                return True
        return False

    @staticmethod
    def _fetch_pdb(pdb_id):
        """ Fetch a PDB with the given id.
        """
        url = f'http://files.rcsb.org/download/{pdb_id}.pdb'
        res = requests.get(url, allow_redirects=True)

        if res.status_code != 200:
            raise PDBFetchError(pdb_id, url)

        return res.content

    def load_pdb(self, file_path):
        """ Loads the receptor file.
        """
        self._pl_complex = PLComplex(file_path)
        self._num_frames += 1

    def load_pdb_id(self, pdb_id):
        """ Download the pdb with given id and save it to a temporary file.
        """
        pdb_str = self._fetch_pdb(pdb_id)
        with tempfile.TemporaryFile() as fp:
            fp.seek(0)
            fp.write(pdb_str)
            self._pl_complex = PLComplex(fp)
        self._num_frames += 1

    def add_frame(self):
        """ Add a new frame to the pharmacophore. """
        self._pharmacophores.append([])
        self._pharmacophores_frames.append(self._num_frames)
        self._num_frames += 1

    def add_points_to_frame(self, point_list, frame):
        """ Add pharmacophoric points from a list to a frame. """
        for point in point_list:
            self._pharmacophores[frame].append(point)

    def add_point(self, point, frame):
        """ Add a pharmacophoric point to a pharmacophore in a specific frame."""
        self._pharmacophores[frame].append(point)

    def remove_point(self, index, frame):
        """ Removes a pharmacophoric point from the pharmacophore at the given frame."""
        self._pharmacophores[frame].pop(index)

    def remove_picked_point(self, view):
        raise NotImplementedError

    def edit_picked_point(self, view):
        raise NotImplementedError

    def add_point_in_picked_location(self, view):
        raise NotImplementedError

    def add_to_view(self, view, frame=0):
        """ Add pharmacophore(s) to a ngl view.
        """
        if isinstance(frame, list):
            raise NotImplementedError
        else:
            for point in self[frame]:
                point.add_to_ngl_view(view)

    def show(self, frame=0):
        """ Shows a 3D representation of the pharmacophore model. """
        view = nv.NGLWidget()
        self.add_to_view(view, frame)
        return view

    def to_json(self, file_name, frame):
        """ Save pharmacophore(s) to a json file.
        """
        data = json_pharmacophoric_elements(self[frame])
        with open(file_name, "w") as fp:
            json.dump(data, fp)

    def to_ligand_scout(self, file_name, frame):
        """ Save a pharmacophore at a given frame to ligand scout format (pml).
        """
        xml_tree = ligandscout_xml_tree(self[frame])
        xml_tree.write(file_name, encoding="UTF-8", xml_declaration=True)

    def to_moe(self, file_name, frame):
        """ Save a pharmacophore at a given frame to moe format (ph4).
        """
        pharmacophore_str = ph4_string(self[frame])
        with open(file_name, "w") as fp:
            fp.write(pharmacophore_str)

    def to_mol2(self, file_name, frame=None):
        """ Save pharmacophore(s) to mol2 file.
        """
        # TODO: save multiple pharmacophores
        if frame is None or isinstance(frame, list):
            raise NotImplementedError
        pharmacophore_data = mol2_file_info([self[frame]])
        with open(file_name, "w") as fp:
            fp.writelines(pharmacophore_data[0])

    def to_rdkit(self, frame):
        """ Transform a pharmacophore at a given frame to a rdkit pharmacophore.
        """
        return rdkit_pharmacophore(self[frame])

    def extract(self, ligand_id, frames=0):
        """ Extract pharmacophore(s) from the receptor. A protein-ligand complex
            can contain multiple ligands or small molecules, pharmacophore(s) is
            extracted only for the selected one.

            Parameters
            ----------
            ligand_id : str
                The id of the ligand whose pharmacophore will be extracted.

            frames : int or list[int] or 'all', optional
                Extract pharmacophores from the given frame(s) of the trajectory.
        """
        if isinstance(frames, int):
            self._pharmacophores.append([])
            self._pharmacophores_frames.append(frames)
            self._num_frames += 1

            ligand_feats = self._pl_complex.ligand_feats_center(ligand_id)
            receptor_feats = self._pl_complex.receptor_feats_center()

            if "hb donor" in ligand_feats:
                self._hb_donor_pharmacophoric_points(
                    ligand_feats["hb donor"], receptor_feats["hb acceptor"], frames)

            if "hb acceptor" in ligand_feats:
                self._hb_acceptor_pharmacophoric_points(
                    ligand_feats["hb acceptor"], receptor_feats["hb donor"], frames)

            if "aromatic ring" in ligand_feats:
                self._aromatic_pharmacophoric_points(
                    ligand_feats["aromatic ring"], receptor_feats["aromatic ring"], frames)

            if "hydrophobicity" in ligand_feats:
                self._hydrophobic_pharmacophoric_points(
                    ligand_feats["hydrophobicity"], receptor_feats["hydrophobicity"], frames)

            if "positive charge" in ligand_feats:
                self._charge_pharmacophoric_points(
                    ligand_feats["positive charge"], receptor_feats["negative charge"],
                    "positive charge", frames)

            if "negative charge" in ligand_feats:
                self._charge_pharmacophoric_points(
                    ligand_feats["negative charge"], receptor_feats["positive charge"],
                    "negative charge", frames)

        else:
            raise NotImplementedError

    def _hb_donor_pharmacophoric_points(self, ligand_centers, receptor_centers, frame):
        """ Compute hydrogen bond donor pharmacophoric points from
            protein-ligand interactions.

            Parameters
            -----------
            ligand_centers : list[puw.Quantity]
                Centroids of the donors in the ligand.

            receptor_centers : list[puw.Quantity]
                Centroids of acceptors in the receptor.

            frame : int
                The frame of the trajectory.

        """
        pass

    def _hb_acceptor_pharmacophoric_points(self, ligand_centers, receptor_centers, frame):
        """ Compute hydrogen bond acceptor pharmacophoric points from
            protein-ligand interactions.

            Parameters
            -----------
            ligand_centers : list[puw.Quantity]
                Centroids of the acceptors in the ligand.

            receptor_centers : list[puw.Quantity]
                Centroids of donors in the receptor.

            frame : int
                The frame of the trajectory.
        """
        pass

    def _aromatic_pharmacophoric_points(self, lig_centers, lig_indices,
                                        rec_centers, rec_indices, frame):
        """ Compute aromatic pharmacophoric points from
            protein-ligand interactions.

            Parameters
            -----------
            lig_centers : list[puw.Quantity]
                Centroids of the aromatic rings in the ligand.

            rec_indices : list[list[int]]
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
                    lig_normal = maths.ring_normal(lig_indices[ii], self._pl_complex.coords, lig_centers[ii])
                    rec_normal = maths.ring_normal(rec_indices[jj], self._pl_complex.coords, rec_centers[jj])
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
                            self._pharmacophores[frame].append(pharma_point)

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
                if maths.points_distance(lig_center, prot_center) < self.HYD_DIST_MAX:
                    centers.append(lig_center)

        points_clustered = self._merge_hydrophobic_points(centers, radius)
        self._pharmacophores[frame] += points_clustered

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
                dist = LigandReceptorPharmacophore._points_distance(centers[ii], centers[jj])
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
                    self._pharmacophores[frame].append(pharma_point)

    @staticmethod
    def _points_distance(coords_1, coords_2):
        """ Returns the distance between two points in 3D space."""
        return np.sqrt(np.sum(np.power(coords_1 - coords_2, 2)))

    def __len__(self):
        return len(self._pharmacophores)

    def __getitem__(self, frame):
        return self._pharmacophores[frame]
