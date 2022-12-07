from openpharmacophore import PharmacophoricPoint
from openpharmacophore.pharmacophore.rdkit_pharmacophore import rdkit_pharmacophore
from openpharmacophore import PLComplex
import openpharmacophore.io as io
from openpharmacophore._private_tools.exceptions import PDBFetchError
from openpharmacophore.utils import maths
import networkx as nx
import numpy as np
import nglview as nv
import pyunitwizard as puw
import json
import re
import requests
import tempfile


class LigandReceptorPharmacophore:
    """ Class to store, and extract pharmacophores from protein-ligand complexes.

        The pharmacophores can be extracted from a pdb file or from a molecular
        dynamics simulation.

    """
    # Values from ligandscout and plip
    HYD_DIST_MAX = puw.quantity(0.5, "nanometers")
    HYD_MERGE_DIST = puw.quantity(0.2, "nanometers")  # value from pharmer

    CHARGE_DIST_MAX = puw.quantity(0.56, "nanometers")

    PISTACK_DIST_MAX = puw.quantity(0.75, "nanometers")
    PISTACK_OFFSET_MAX = puw.quantity(0.20, "nanometers")
    PISTACK_ANG_DEV = 30  # degrees

    def __init__(self):
        self._pharmacophores = []
        self._pharmacophores_frames = []  # Contains the frame to which each pharmacophore belongs
        self._num_frames = 0

        self._pl_complex = None

    @property
    def num_frames(self):
        return self._num_frames

    @property
    def receptor(self):
        """ Return the receptor.

            Returns
            -------
            PLComplex
        """
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

    def load_receptor(self, file_path):
        """ Loads the receptor (or protein-ligand complex) file.

            Can be a pdb or a trajectory file format.

            Parameters
            ----------
            file_path : str

        """
        self._pl_complex = PLComplex(file_path)

    def load_pdb_id(self, pdb_id):
        """ Download the pdb with given id and save it to a temporary file.
        """
        pdb_str = self._fetch_pdb(pdb_id)
        with tempfile.TemporaryFile() as fp:
            fp.seek(0)
            fp.write(pdb_str)
            self._pl_complex = PLComplex(fp)

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
        """ Add pharmacophoric points to a ngl view.

            Parameters
            ----------
            view : nglview.NGLWidget
                A view object where the pharmacophoric points will be added.
            frame : int
                Adds the points of this frame.
        """
        for point in self[frame]:
            point.add_to_ngl_view(view)

    def show(self, frame=0, ligand=True, receptor=True,
             points=True, indices=None):
        """ Shows a 3D representation of the pharmacophore model.

            Parameters
            ----------
            frame : int
                Which frame to show
            ligand : bool, optional
                Whether to show the ligand.
            receptor : bool, optional
                Whether to show the receptor.
            points : bool, optional
                Whether to show pharmacophoric points.
            indices : np.ndarray or list[int]
                A list of the indices of the atoms that will be shown.
        """
        view = nv.NGLWidget()
        if not any([ligand, receptor, points]):
            return view

        if ligand and frame != 0:
            conformer = self.receptor.get_lig_conformer(frame)
        else:
            conformer = self.receptor.ligand

        if indices is None:

            if receptor and ligand:
                view = nv.show_mdtraj(self.receptor.traj.atom_slice(
                    self.receptor.receptor_indices)[frame])
                view.add_component(conformer)
            elif receptor:
                view = nv.show_mdtraj(self.receptor.traj.atom_slice(
                    self.receptor.receptor_indices))
            elif ligand:
                view = nv.show_rdkit(conformer)
        else:
            traj = self.receptor.slice_traj(indices, frame)
            view = nv.show_mdtraj(traj)
            view.representations = [{
                "type": "ball+stick",
                "params": {
                    "sele": "all",
                }
            }]
            if ligand:
                view.add_component(conformer)

        if points:
            self.add_to_view(view, frame)
        return view

    def to_json(self, file_name, frame):
        """ Save pharmacophore(s) to a json file.
        """
        data = io.json_pharmacophoric_elements(self[frame])
        with open(file_name, "w") as fp:
            json.dump(data, fp)

    def to_ligand_scout(self, file_name, frame):
        """ Save a pharmacophore at a given frame to ligand scout format (pml).
        """
        xml_tree = io.ligandscout_xml_tree(self[frame])
        xml_tree.write(file_name, encoding="UTF-8", xml_declaration=True)

    def to_moe(self, file_name, frame):
        """ Save a pharmacophore at a given frame to moe format (ph4).
        """
        pharmacophore_str = io.ph4_string(self[frame])
        with open(file_name, "w") as fp:
            fp.write(pharmacophore_str)

    def to_mol2(self, file_name, frame=None):
        """ Save pharmacophore(s) to mol2 file.
        """
        # TODO: save multiple pharmacophores
        if frame is None or isinstance(frame, list):
            raise NotImplementedError
        pharmacophore_data = io.mol2_file_info([self[frame]])
        with open(file_name, "w") as fp:
            fp.writelines(pharmacophore_data[0])

    def to_rdkit(self, frame):
        """ Transform a pharmacophore at a given frame to a rdkit pharmacophore.
        """
        return rdkit_pharmacophore(self[frame])

    def extract(self, ligand_id, frames=0, smiles="", features=None,
                add_hydrogens=True):
        """ Extract pharmacophore(s) from the receptor. A protein-ligand complex
            can contain multiple ligands or small molecules, pharmacophore(s) is
            extracted only for the selected one.

            Parameters
            ----------
            ligand_id : str
                The id of the ligand whose pharmacophore will be extracted.

            frames : int or list[int] or 'all', optional
                Extract pharmacophores from the given frame(s) of the trajectory.

            smiles : str, optional
                The smiles of the ligand.

            features : list[str], optional
                A list of the chemical features that will be used in the extraction.

            add_hydrogens : bool
                Whether to add hydrogens to the ligand and the receptor. Necessary
                to extrac hydrogen bond pharmacophoric points.
        """
        if isinstance(frames, int):
            frames = [frames]

        if features is None:
            features = PharmacophoricPoint.get_valid_features()

        pl: PLComplex
        pl = self._pl_complex
        pl.prepare(lig_id=ligand_id, smiles=smiles, add_hydrogens=add_hydrogens)

        for frame in frames:
            self.add_frame()

            if "hydrophobicity" in features:
                lig_hyd_centers, _ = pl.ligand_features("hydrophobicity", frame)
                if len(lig_hyd_centers) > 0:
                    rec_hyd_centers, _ = pl.receptor_features("hydrophobicity", frame)
                    self._hydrophobic_pharmacophoric_points(
                        lig_hyd_centers, rec_hyd_centers, frame
                    )

            if "positive charge" in features:
                lig_pos_charge_cent, _ = pl.ligand_features("positive charge", frame)
                if len(lig_pos_charge_cent) > 0:
                    rec_neg_charge_cent, _ = pl.receptor_features("negative charge", frame)
                    self._charge_pharmacophoric_points(
                        lig_pos_charge_cent, rec_neg_charge_cent,
                        "positive charge", frame
                    )

            if "negative charge" in features:
                lig_neg_charge_cent, _ = pl.ligand_features("negative charge", frame)
                if len(lig_neg_charge_cent) > 0:
                    rec_pos_charge_cent, _ = pl.receptor_features("positive charge", frame)
                    self._charge_pharmacophoric_points(
                        lig_neg_charge_cent, rec_pos_charge_cent,
                        "negative charge", frame
                    )

            if "aromatic ring" in features:
                lig_aro_cent, lig_aro_ind = pl.ligand_features("aromatic ring", frame)
                if len(lig_aro_cent) > 0:
                    rec_aro_cent, rec_aro_ind = pl.receptor_features("aromatic ring", frame)
                    self._aromatic_pharmacophoric_points(
                        lig_aro_cent, lig_aro_ind,
                        rec_aro_cent, rec_aro_ind, frame
                    )

            h_bonds = None
            if "hb donor" in features:
                h_bonds = pl.hbond_indices(frame)
                self._hbond_donor_pharmacophoric_points(h_bonds, frame)

            if "hb acceptor" in features:
                if h_bonds is None:
                    h_bonds = pl.hbond_indices(frame)
                self._hbond_acceptor_pharmacophoric_points(h_bonds, frame)

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
                self._pharmacophores[frame].append(pharma_point)

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
                self._pharmacophores[frame].append(pharma_point)

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
                dist = maths.points_distance(lig_center, prot_center)
                if dist < self.HYD_DIST_MAX:
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
                    self._pharmacophores[frame].append(pharma_point)

    def __len__(self):
        return len(self._pharmacophores)

    def __getitem__(self, frame):
        return self._pharmacophores[frame]
