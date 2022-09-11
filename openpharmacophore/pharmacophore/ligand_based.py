from .pharmacophoric_point import distance_between_pharmacophoric_points, PharmacophoricPoint
from .pharmacophore import Pharmacophore
from .rdkit_pharmacophore import rdkit_pharmacophore
from ..io import (json_pharmacophoric_elements, ligandscout_xml_tree,
                  mol2_file_info, ph4_string)
from ..io import (load_json_pharmacophore, load_mol2_pharmacophoric_points,
                  pharmacophoric_points_from_ph4_file, read_ligandscout)
from .._private_tools.exceptions import InvalidFileFormat
import numpy as np
import nglview as nv
import pyunitwizard as puw
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
        """
               Load a pharmacophore from a file.

               Parameters
               ---------
               file_name : str
                   Name of the file containing the pharmacophore

       """
        if file_name.endswith(".json"):
            return load_json_pharmacophore(file_name)[0]
        elif file_name.endswith(".mol2"):
            return load_mol2_pharmacophoric_points(file_name)[0]
        elif file_name.endswith(".pml"):
            return read_ligandscout(file_name)
        elif file_name.endswith("ph4"):
            return pharmacophoric_points_from_ph4_file(file_name)
        else:
            raise InvalidFileFormat(file_name.split(".")[-1])

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

    @staticmethod
    def get_picked_point_index(view):
        """ Get the index of a point picked in a view in the
            pharmacophore.

            Parameters
            ----------
            view : nglview.NGLWidget
                View where the pharmacophore will be added.

            Returns
            -------
            index: int or None
                The index or None if something invalid was selected.
        """
        if len(view.picked) != 1:
            # An atom or nothing was selected
            return

        picked_index = view.picked["component"]
        # We'll assume that the view components are the molecule(s) followed by
        # the pharmacophoric points
        points_start = None
        for ii in range(len(view._ngl_component_names)):
            if view._ngl_component_names[ii] == "nglview.shape.Shape":
                points_start = ii
                break

        if points_start is None:
            return

        return picked_index - points_start

    def remove_picked_point(self, view):
        """ Remove a pharmacophoric point selected in a ngl view.

            If nothing is selected or the selected thing is not a
            pharmacophoric point, it doesn't change the view.

           Parameters
           ----------
           view : nglview.NGLWidget
               View where the pharmacophore will be added.
        """
        index = self.get_picked_point_index(view)
        if index is not None and index < len(self):
            self.remove_point(index)

    def edit_picked_point(self, view, center, radius):
        """ Remove a pharmacophoric point selected in a ngl view.

           If nothing is selected or the selected thing is not a
           pharmacophoric point, it doesn't change the view.

          Parameters
          ----------
          view : nglview.NGLWidget
              View where the pharmacophore will be added.

          center: Quantity
              The new center of the pharmacophoric point

          radius: Quantity:
              The new radius of the pharmacophric point
        """
        index = self.get_picked_point_index(view)
        if index is not None and index < len(self):
            self[index].center = center
            self[index].radius = radius

    def add_point_in_picked_location(self, view, feat_name, radius):
        """ Add a pharmacophoric point in the selected atom of a ngl view.

            Parameters
              ----------
              view : nglview.NGLWidget
                  View where the pharmacophore will be added.

              feat_name: str
                  The feature type of the pharmacophoric point

              radius: Quantity:
                  The new radius of the pharmacophric point

        """
        try:
            atom = view.picked["atom1"]
        except KeyError:
            return

        center = puw.quantity([atom["x"], atom["y"], atom["z"]], "angstroms")
        new_point = PharmacophoricPoint(feat_name, center, radius)
        self._points.append(new_point)

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
        """ Save the pharmacophore to a ligand scout file (pml).

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        xml_tree = ligandscout_xml_tree(self.pharmacophoric_points)
        xml_tree.write(file_name, encoding="UTF-8", xml_declaration=True)

    def to_moe(self, file_name):
        """ Save the pharmacophore to a moe file (ph4).

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        pharmacophore_str = ph4_string(self.pharmacophoric_points)
        with open(file_name, "w") as fp:
            fp.write(pharmacophore_str)

    def to_mol2(self, file_name):
        """ Save the pharmacophore to a mol2 file.

            Parameters
            ----------
            file_name: str
                Name of the json file.
        """
        pharmacophore_data = mol2_file_info([self.pharmacophoric_points])
        with open(file_name, "w") as fp:
            fp.writelines(pharmacophore_data[0])

    def to_rdkit(self):
        """ Returns a rdkit pharmacophore with the pharmacophoric_points from the original pharmacophore.

            rdkit pharmacophores do not store the pharmacophoric_points radii, so they are returned as well.

            Returns
            -------
            rdkit_pharmacophore : rdkit.Chem.Pharm3D.Pharmacophore
                The rdkit pharmacophore.

            radii : list[float]
                List with the radius in angstroms of each pharmacophoric point.
        """
        return rdkit_pharmacophore(self.pharmacophoric_points)

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
