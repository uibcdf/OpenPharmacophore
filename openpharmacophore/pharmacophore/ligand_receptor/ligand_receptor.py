from openpharmacophore import PharmacophoricPoint, Pharmacophore
from openpharmacophore.molecular_systems import ComplexBindingSite, Ligand
from openpharmacophore.utils import maths
from openpharmacophore.constants import FEAT_TYPES
import networkx as nx
import numpy as np
import pyunitwizard as puw


class LigandReceptorPharmacophore:
    """ Class to store, and extract pharmacophores from protein-ligand complexes.

        The pharmacophores can be extracted from a PDB structure or from a molecular
        dynamics trajectory.

        Parameters
        ----------
        binding_site : ComplexBindingSite
            Binding site from which the pharmacophore will be extracted

        ligand : Ligand
            The ligand in the binding site.

    """
    # Values from ligandscout and plip
    HYD_DIST_MAX = puw.quantity(0.5, "nanometers")
    HYD_MERGE_DIST = puw.quantity(0.2, "nanometers")  # value from pharmer

    CHARGE_DIST_MAX = puw.quantity(0.56, "nanometers")

    PISTACK_DIST_MAX = puw.quantity(0.75, "nanometers")
    PISTACK_OFFSET_MAX = puw.quantity(0.20, "nanometers")
    PISTACK_ANG_DEV = 30  # degrees

    HBOND_DIST_MAX = puw.quantity(0.25, "nanometers")
    HBOND_ANG_MIN = puw.quantity(120, "degree")  # degrees

    def __init__(self, binding_site, ligand):
        self._pharmacophores = []  # type: list[Pharmacophore]
        self._bsite = binding_site
        self._ligand = ligand

    def _hydrophobic_pharmacophoric_points(self, ligand_feats, receptor_feats):
        """ Get hydrophobic pharmacophoric points.

            Parameters
            ----------
            ligand_feats : ChemFeatContainer
            receptor_feats : ChemFeatContainer

            Returns
            -------
            list [PharmacophoricPoint]

        """
        hyd_points = self._points_from_distance_rule(
            ligand_feats.hydrophobic, receptor_feats.hydrophobic, "hydrophobicity", self.HYD_DIST_MAX
        )
        return self._merge_hydrophobic_points(hyd_points, radius=puw.quantity(1.0, "angstroms"))

    def extract(self, frames=None, feat_types=FEAT_TYPES):
        """ Extract one or multiple pharmaophores at different frames

            Parameters
            ----------
            frames : Iterable[int], optional
                Iterable with the frames for which pharmacophores will be extracted.

            feat_types : set[str], optional
                Only pharmacophoric points with this feature types will be extracted.
        """
        if frames is None:
            frames = range(1)

        for fr in frames:
            self._extract_single_frame(fr, feat_types)

    def _extract_single_frame(self, frame, feat_types=FEAT_TYPES):
        """ Extract a pharmacophore at a given frame.

            Parameters
            ----------
            frame : int
                The frame of the trajectory for which the pharmacophore will be
                extracted.

            feat_types : set[str], optional
                Only pharmacophoric points with this feature types will be extracted.
        """
        ligand_feats = self._ligand.get_chem_feats_with_directionality(frame)
        receptor_feats = self._bsite.get_chem_feats(frame)

        points = []
        if "positive charge" in feat_types:
            points += self._points_from_distance_rule(
                ligand_feats.positive, receptor_feats.negative, "positive charge", self.CHARGE_DIST_MAX
            )
        if "negative charge" in feat_types:
            points += self._points_from_distance_rule(
                ligand_feats.negative, receptor_feats.positive, "negative charge", self.CHARGE_DIST_MAX
            )
        if "aromatic ring" in feat_types:
            points += self._aromatic_pharmacophoric_points(ligand_feats.aromatic, receptor_feats.aromatic)
        if "hb donor" in feat_types:
            points += self._hbond_pharmacophoric_points(
                ligand_feats.donor, receptor_feats.acceptor, donors_in_ligand=True
            )
        if "hb acceptor" in feat_types:
            points += self._hbond_pharmacophoric_points(
                receptor_feats.donor, ligand_feats.acceptor, donors_in_ligand=False
            )
        if "hydrophobicity" in feat_types:
            points += self._hydrophobic_pharmacophoric_points(ligand_feats, receptor_feats)

        pharma = Pharmacophore(points, ref_struct=frame, ref_mol=0)
        self._pharmacophores.append(pharma)

    @staticmethod
    def _hbond_angle(don_acc_dist, don_hyd_dist, acc_hyd_dist):
        """  Compute hydrogen bonding angle.

            Parameters
            ----------
            don_acc_dist : QuantityLike
                Distance from donor to acceptor

            don_hyd_dist : QuantityLike
                Distance from donor to hydrogen

            acc_hyd_dist : QuantityLike
               Distance from acceptor to hydrogen

            Returns
            -------
            QuantityLike
        """
        # Calculate angle using law of cosines
        cos = (don_hyd_dist**2 + acc_hyd_dist**2 - don_acc_dist**2) / (2 * don_hyd_dist * acc_hyd_dist)
        return np.degrees(np.arccos(cos))

    @staticmethod
    def _hbond_pharmacophoric_points(donors, acceptors, donors_in_ligand):
        """ Compute hydrogen bonds pharmacophoric points centroids and
            directions vectors.

            Parameters
            ----------
            donors : list[HBDonor]
                List with hydrogen bond donors.

            acceptors : list[ChemFeat]
                List with hydrogen bond acceptors.

            donors_in_ligand : bool
                Whether the donors are part of the ligand or of the receptor

            Returns
            -------
            points : list[PharmacophoricPoint]

        """
        radius = puw.quantity(1.0, "angstroms")

        points = []
        for don in donors:
            for acc in acceptors:
                # Distance between H-Acceptor
                acc_hyd_dist = maths.points_distance(don.hyd, acc.coords)
                if acc_hyd_dist < LigandReceptorPharmacophore.HBOND_DIST_MAX:
                    don_acc_dist = maths.points_distance(don.coords, acc.coords)
                    don_hyd_dist = maths.points_distance(don.coords, don.hyd)
                    angle = LigandReceptorPharmacophore._hbond_angle(
                        don_acc_dist, don_hyd_dist, acc_hyd_dist
                    )
                    if angle > LigandReceptorPharmacophore.HBOND_ANG_MIN:
                        direction = puw.get_value(acc.coords - don.hyd)
                        if donors_in_ligand:
                            points.append(
                                PharmacophoricPoint(
                                    "hb donor", don.coords, radius, direction=direction)
                            )
                        else:
                            points.append(
                                PharmacophoricPoint(
                                    "hb acceptor", acc.coords, radius, direction=direction)
                            )
        return points

    @staticmethod
    def _aromatic_angle_exceeds_deviation(angle):
        ang_dev = LigandReceptorPharmacophore.PISTACK_ANG_DEV
        return not (0 <= angle <= ang_dev or 90 - ang_dev <= angle <= 90 + ang_dev)

    @staticmethod
    def _aromatic_pharmacophoric_points(ligand_feats, receptor_feats):
        """ Compute aromatic pharmacophoric points from
            protein-ligand interactions.

            Parameters
            -----------
            ligand_feats : list[ChemFeat]
                Aromatic chemical features of the ligand

            receptor_feats : list[ChemFeat]
                Aromatic chemical features of the receptor

            Returns
            -------
            points : list[PharmacophoricPoint]
                List of aromatic pharmacophoric points

        """
        # Calculate pistack interactions and create pharmacophoric points
        radius = puw.quantity(1.0, "angstroms")
        points = []
        for ii in range(len(ligand_feats)):
            for jj in range(len(receptor_feats)):
                lig_center = ligand_feats[ii].coords
                rec_center = receptor_feats[jj].coords
                dist = maths.points_distance(rec_center, lig_center)
                if dist <= LigandReceptorPharmacophore.PISTACK_DIST_MAX:
                    # Calculate deviation from ideal angle by taking the angle between the normals
                    # defined by the planes of each ring
                    lig_normal = ligand_feats[ii].normal
                    rec_normal = receptor_feats[jj].normal

                    angle = maths.angle_between_normals(lig_normal, rec_normal)
                    assert 0 <= angle <= 360, f"Angle is {angle}"

                    if not LigandReceptorPharmacophore._aromatic_angle_exceeds_deviation(angle):

                        # Project ring centers into the other plane and calculate offset
                        rec_proj = maths.point_projection(
                            lig_normal, lig_center, rec_center
                        )
                        lig_proj = maths.point_projection(
                            rec_normal, rec_center, lig_center
                        )
                        offset = min(maths.points_distance(lig_proj, rec_center),
                                     maths.points_distance(rec_proj, lig_center))

                        if offset <= LigandReceptorPharmacophore.PISTACK_OFFSET_MAX:
                            direction = puw.get_value(rec_center - lig_center)
                            points.append(PharmacophoricPoint(
                                "aromatic ring", lig_center, radius, direction
                            ))
        return points

    @staticmethod
    def _merge_hydrophobic_points(points, radius):
        """ Merge group of hydrophobic points close to each other.

            Parameters
            ----------
            points : list[PharmacophoricPoint]
            radius : puw.Quantity

            Returns
            -------
            list[PharmacophoricPoint]
        """
        # Create a graph of the hydrophobic features, were each node represents
        # a feature and an edge is added between two nodes only if their distance
        # is < HYD_MERGE_DIST.
        hyd_graph = nx.Graph()
        for ii in range(len(points)):
            hyd_graph.add_node(ii)

        for ii in range(len(points)):
            for jj in range(ii + 1, len(points)):
                dist = maths.points_distance(points[ii].center, points[jj].center)
                if dist <= LigandReceptorPharmacophore.HYD_MERGE_DIST:
                    hyd_graph.add_edge(ii, jj)

        # Find each maximum clique and group all nodes within a clique
        cliques_iter = nx.find_cliques(hyd_graph)
        merged = []
        for clique in cliques_iter:
            clique_centers = [points[ii].center for ii in clique]
            center = np.mean(np.stack(clique_centers), axis=0)
            pharma_point = PharmacophoricPoint("hydrophobicity", center, radius)
            merged.append(pharma_point)
        return merged

    @staticmethod
    def _points_from_distance_rule(
            ligand_feats, receptor_feats, feat_type, max_dist
    ):
        """ Given a set of ligand and receptor chem feats, creates pharmacophoric
            points if the ligand and receptor feats are below the maximum distance.

            This method can be used to obtain hydrophobic, positive charge and
            negative charge pharmacophoric points.

        Parameters
        ----------
        ligand_feats : list[ChemFeat]
            Chemical features of the ligand

        receptor_feats : list[ChemFeat]
            Chemical features of the receptor

        feat_type : str
            Name of the feature type of the pharmacophoric points.

        max_dist : QuantityLike
            Maximum distance between a receptor and a ligand feature.

        Returns
        -------
        points : list[PharmacophoricPoint]

        """
        radius = puw.quantity(1.0, "angstroms")
        points = []
        for lig_feat in ligand_feats:
            for rec_feat in receptor_feats:
                dist = maths.points_distance(lig_feat.coords, rec_feat.coords)
                if dist < max_dist:
                    pharma_point = PharmacophoricPoint(feat_type, lig_feat.coords, radius)
                    points.append(pharma_point)
        return points

    def __len__(self):
        return len(self._pharmacophores)

    def __getitem__(self, frame):
        """
            Returns
            -------
            Pharmacophore
        """
        return self._pharmacophores[frame]
