from .pharmacophoric_point import distance_between_pharmacophoric_points
from .pharmacophore import Pharmacophore
from ..io import (json_pharmacophoric_elements, ligandscout_xml_tree,
                  mol2_file_info, ph4_string)
from ..io import (load_json_pharmacophore, load_mol2_pharmacophoric_points,
                  pharmacophoric_points_from_ph4_file, read_ligandscout)
import numpy as np
import nglview as nv
import pyunitwizard as puw
from rdkit import Geometry
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm3D import Pharmacophore as rdkitPharmacophore
import json


class LigandBasedPharmacophore(Pharmacophore):
    """ Class to store and extract pharmacophores from a set of ligands.

        Parameters
        ----------
        ligands : str or list[rdkit.Mol]
            A file with ligands or a list of molecules
    """

    def __init__(self, ligands):
        if isinstance(ligands, str):
            self._points = self.from_file(ligands)
        else:
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

    def from_file(self, file_name):

        if file_name.endswith(".json"):
            return load_json_pharmacophore(file_name)[0]
        elif file_name.endswith(".mol2"):
            return load_mol2_pharmacophoric_points(file_name)[0]
        elif file_name.endswith(".pml"):
            return read_ligandscout(file_name)
        elif file_name.endswith("ph4"):
            return pharmacophoric_points_from_ph4_file(file_name)
        else:
            raise ValueError

    def add_point(self, point):
        """ Adds a pharmacophoric point.

            Parameters
            ----------
            point : PharmacophoricPoint
                The pharmacophoric point that will be added.
        """
        self._points.append(point)

    def remove_point(self, index):
        """ Removes a pharmacophoric point from the pharmacophore.

            Parameters
            ----------
            index: int
                The index of the pharmacophoric point.
        """
        self._points.pop(index)

    def remove_picked_point(self, *args, **kwargs):
        pass

    def edit_picked_point(self, *args, **kwargs):
        pass

    def add_point_in_picked_location(self, *args, **kwargs):
        pass

    def add_to_view(self, view, palette=None, opacity=0.5):
        """Add the pharmacophore representation to a view from NGLView.

           Parameters
           ----------
           view : nglview.NGLWidget
               View where the pharmacophore will be added.

           palette : dict[str, str], optional
               Dictionary with a color for each feature type.

           opacity : float
                The level of opacity of the points. Must be a number between 0 and 1.

        """
        for point in self.pharmacophoric_points:
            point.add_to_ngl_view(view, palette, opacity)

    def show(self, palette=None):
        """ Show the pharmacophore model.

        Parameters
        ----------
        palette : str or dict.
            Color palette name or dictionary. (Default: 'openpharmacophore')

        Returns
        -------
        nglview.NGLWidget
            A nglview.NGLWidget with the 'view' of the pharmacophoric model.
        """
        view = nv.NGLWidget()
        self.add_to_view(view, palette=palette)

        return view

    def to_json(self, file_name):
        """ Save the pharmacophore to a json file.

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        data = json_pharmacophoric_elements(self.pharmacophoric_points)
        with open(file_name, "w") as fp:
            json.dump(data, fp)

    def to_ligand_scout(self, file_name):
        xml_tree = ligandscout_xml_tree(self.pharmacophoric_points)
        xml_tree.write(file_name, encoding="UTF-8", xml_declaration=True)

    def to_moe(self, file_name):
        pharmacophore_str = ph4_string(self.pharmacophoric_points)
        with open(file_name, "w") as fp:
            fp.write(pharmacophore_str)

    def to_mol2(self, file_name):
        pharmacophore_data = mol2_file_info(self)
        with open(file_name, "w") as fp:
            fp.writelines(pharmacophore_data[0])

    def to_rdkit(self):
        """ Returns a rdkit pharmacophore with the pharmacophoric_points from the original pharmacophore.

            rdkit pharmacophores do not store the pharmacophoric_points radii, so they are returned as well.

            Returns
            -------
            rdkit_pharmacophore : rdkit.Chem.Pharm3D.Pharmacophore
                The rdkit pharmacophore.

            radii : list of float
                List with the radius in angstroms of each pharmacophoric point.
        """
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
        """ Compute the distance matrix of the pharmacophore.

            Returns
            -------
            dis_matrix : np.ndarray of shape(num_points, num_points)
                The distance matrix.
        """
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
        """ Extracts a pharmacophore from a set of ligands.
        """
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
