from .pharmacophoric_point import distance_between_pharmacophoric_points
from openpharmacophore.pharmacophore.pharmacophore import Pharmacophore
import numpy as np
import pyunitwizard as puw
from rdkit import Geometry
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm3D import Pharmacophore as rdkitPharmacophore


class LigandBasedPharmacophore(Pharmacophore):

    def __init__(self, ligands):
        self._points = self._extract_pharmacophore(ligands)

    @property
    def pharmacophoric_points(self):
        return self._points

    @pharmacophoric_points.setter
    def pharmacophoric_points(self, points):
        self._points = points

    @property
    def num_points(self, *args, **kwargs):
        return len(self._points)

    @classmethod
    def from_file(cls, file_name):
        pass

    def add_point(self, *args, **kwargs):
        pass

    def remove_point(self, *args, **kwargs):
        pass

    def remove_picked_point(self, *args, **kwargs):
        pass

    def edit_picked_point(self, *args, **kwargs):
        pass

    def add_point_in_picked_location(self, *args, **kwargs):
        pass

    def to_json(self, *args, **kwargs):
        pass

    def add_to_view(self, *args, **kwargs):
        pass

    def show(self, *args, **kwargs):
        pass

    def to_ligand_scout(self, *args, **kwargs):
        pass

    def to_moe(self, file_name):
        pass

    def to_pharmagist(self, file_name):
        pass

    def to_rdkit(self):
        rdkit_element_name = {
            "aromatic ring": "Aromatic",
            "hydrophobicity": "Hydrophobe",
            "hb acceptor": "Acceptor",
            "hb donor": "Donor",
            "positive charge": "PosIonizable",
            "negative charge": "NegIonizable",
        }

        points = []
        radii = []

        for element in self:
            feat_name = rdkit_element_name[element.feature_name]
            center = puw.get_value(element.center, to_unit="angstroms")
            center = Geometry.Point3D(center[0], center[1], center[2])
            points.append(ChemicalFeatures.FreeChemicalFeature(
                feat_name,
                center
            ))
            radius = puw.get_value(element.radius, to_unit="angstroms")
            radii.append(radius)

        rdkit_pharmacophore = rdkitPharmacophore.Pharmacophore(points)
        return rdkit_pharmacophore, radii

    def distance_matrix(self):
        n_pharmacophoric_points = self.num_points
        dis_matrix = np.zeros((n_pharmacophoric_points, n_pharmacophoric_points))

        for ii in range(n_pharmacophoric_points):
            for jj in range(ii, n_pharmacophoric_points):
                if ii == jj:
                    dis_matrix[ii, jj] = 0
                else:
                    distance = distance_between_pharmacophoric_points(
                        self[ii], self[jj])
                    dis_matrix[ii, jj] = distance
                    dis_matrix[jj, ii] = distance

        return dis_matrix

    @staticmethod
    def _extract_pharmacophore(ligands):
        return []

    def __len__(self):
        return len(self._points)

    def __getitem__(self, index):
        return self._points[index]

    def __eq__(self, other):
        """ Check equality between pharmacophores.

            Assumes that pharmacophoric points are sorted equally in both pharmacophores.
        """
        if isinstance(other, type(self)) and self.num_points == other.num_points:
            for ii in range(self.num_points):
                if not self.pharmacophoric_points[ii] == other.pharmacophoric_points[ii]:
                    return False
            return True
        return False

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n_pharmacophoric_points: {self.num_points})"
