# OpenPharmacophore
from openpharmacophore._private_tools.exceptions import InvalidFeatureError, InvalidFileFormat
from openpharmacophore.io import (from_pharmer, from_moe, from_ligandscout, _ligandscout_xml_tree, _moe_ph4_string,
                                  to_pharmagist, _pharmer_dict)
from openpharmacophore import PharmacophoricPoint
from openpharmacophore.utils.discretize import discretize
from openpharmacophore.utils.bisection import insort_right
from openpharmacophore.pharmacophore.pharmacophoric_point import distance_between_pharmacophoric_points
from openpharmacophore.pharmacophore.color_palettes import get_color_from_palette_for_feature
# Third party
import networkx as nx
import nglview as nv
import numpy as np
import pyunitwizard as puw
from rdkit import Geometry, RDLogger
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm3D import Pharmacophore as rdkitPharmacophore
# Standard library
import json
from typing import List, Dict, Tuple

RDLogger.DisableLog('rdApp.*')  # Disable rdkit warnings


class Pharmacophore:
    """ Native object for pharmacophores.

    Parameters
    ----------

    pharmacophoric_points : list openpharmacophore.PharmacophoricPoint
        List of pharmacophoric points. 

    Attributes
    ----------
    _pharmacophoric_points : list of openpharmacophore.PharmacophoricPoint
        Private attribute, should not be modified directly. To add or remove pharmacophoric 
        points use the respective methods

    n_pharmacophoric_points : int
        Number of pharmacophoric points

    """

    def __init__(self, pharmacophoric_points: List[PharmacophoricPoint]) -> None:

        if pharmacophoric_points is not None:
            self._pharmacophoric_points = sorted(pharmacophoric_points, key=lambda p: p.short_name)
        self.n_pharmacophoric_points = len(pharmacophoric_points)

    @classmethod
    def from_file(cls, file_name: str) -> "Pharmacophore":
        """
        Class method to load a pharmacophore from a file.

        Parameters
        ---------
        file_name : str
            Name of the file containing the pharmacophore

        """
        file_extension = file_name.split(".")[-1]
        if file_extension == "json":
            points, _, _ = from_pharmer(file_name, False)

        elif file_extension == "ph4":
            points = from_moe(file_name)

        elif file_extension == "pml":
            points = from_ligandscout(file_name)

        else:
            raise InvalidFileFormat(f"Invalid file format, \"{file_name}\" is not a supported file format")

        return cls(pharmacophoric_points=points)

    def add_to_NGLView(self, view: nv.NGLWidget, palette: str = 'openpharmacophore') -> None:
        """Add the pharmacophore representation to a view (NGLWidget) from NGLView.

        Each pharmacophoric element is added to the NGLWidget as a new component.

        Parameters
        ----------
        view : nglview.NGLWidget
            View as NGLView widget where the pharmacophore will be added.

        palette : str or dict
            Color palette name or dictionary. (Default: 'openpharmacophore')

        """
        first_element_index = len(view._ngl_component_ids)
        for ii, element in enumerate(self._pharmacophoric_points):
            # Add Spheres
            center = puw.get_value(element.center, to_unit="angstroms").tolist()
            radius = puw.get_value(element.radius, to_unit="angstroms")
            feature_color = get_color_from_palette_for_feature(element.feature_name, color_palette=palette)
            label = f"{element.feature_name}_{ii}"
            view.shape.add_sphere(center, feature_color, radius, label)

            # Add vectors
            if element.has_direction:
                label = f"{element.feature_name}_vector"
                if element.feature_name == "hb acceptor":
                    end_arrow = puw.get_value(
                        element.center - 2 * radius * puw.quantity(element.direction, "angstroms"),
                        to_unit='angstroms').tolist()
                    view.shape.add_arrow(end_arrow, center, feature_color, 0.2, label)
                else:
                    end_arrow = puw.get_value(
                        element.center + 2 * radius * puw.quantity(element.direction, "angstroms"),
                        to_unit='angstroms').tolist()
                    view.shape.add_arrow(center, end_arrow, feature_color, 0.2, label)

        # Add opacity to spheres
        last_element_index = len(view._ngl_component_ids)
        for jj in range(first_element_index, last_element_index):
            view.update_representation(component=jj, opacity=0.8)

    def show(self, palette: str = 'openpharmacophore'):
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
        self.add_to_NGLView(view, palette=palette)

        return view

    def add_point(self, pharmacophoric_point: PharmacophoricPoint) -> None:
        """ Adds a pharmacophoric point.

            Parameters
            ----------
            pharmacophoric_point : PharmacophoricPoint
                The pharmacophoric point that will be added.
        """
        insort_right(self._pharmacophoric_points, pharmacophoric_point, key=lambda p: p.short_name)
        self.n_pharmacophoric_points += 1

    @property
    def pharmacophoric_points(self) -> List[PharmacophoricPoint]:
        return self._pharmacophoric_points

    def remove_point(self, index: int) -> None:
        """Remove a pharmacophoric point from the pharmacophore.
        """
        self._pharmacophoric_points.pop(index)
        self.n_pharmacophoric_points -= 1

    def remove_points(self, indices: List[int]) -> None:
        """ Remove a list of pharmacophoric points from the pharmacophore.
        """
        for index in indices:
            if index < 0 or index > self.n_pharmacophoric_points - 1:
                raise IndexError

        new_pharmacophoric_points = []
        for ii, point in enumerate(self._pharmacophoric_points):
            if ii not in indices:
                new_pharmacophoric_points.append(point)

        self._pharmacophoric_points = new_pharmacophoric_points
        self.n_pharmacophoric_points = len(self._pharmacophoric_points)

    def remove_feature(self, feat_type: str) -> None:
        """ Remove an especific feature type from the pharmacophore pharmacophoric_points list

            Parameters
            ----------
            feat_type : str
                Name or type of the feature to be removed.
        """
        feats = PharmacophoricPoint.get_valid_features()
        if feat_type not in feats:
            raise InvalidFeatureError(f"Cannot remove feature. \"{feat_type}\" is not a valid feature type")

        temp_pharmacophoric_points = [element for element in self._pharmacophoric_points if
                                      element.feature_name != feat_type]
        self._pharmacophoric_points = temp_pharmacophoric_points
        self.n_pharmacophoric_points = len(self._pharmacophoric_points)

    def _reset(self) -> None:
        """Private method to reset all attributes to default values.
        """
        self._pharmacophoric_points.clear()
        self.n_pharmacophoric_points = 0

    def to_ligandscout(self, file_name: str) -> None:
        """Method to export the pharmacophore to the ligandscout compatible format.

            Parameters
            ----------
            file_name : str
                Name of file to be written with the ligandscout format of the pharmacophore.

        """
        tree = _ligandscout_xml_tree(self._pharmacophoric_points)
        tree.write(file_name, encoding="UTF-8", xml_declaration=True)

    def to_pharmer(self, file_name: str) -> None:
        """Method to export the pharmacophore to the pharmer compatible format.

            Parameters
            ----------
            file_name : str
                Name of file to be written with the pharmer format of the pharmacophore.
        """
        pharmacophore_dict = _pharmer_dict(self.pharmacophoric_points)

        with open(file_name, "w") as outfile:
            json.dump(pharmacophore_dict, outfile)

    def to_pharmagist(self, file_name: str) -> None:
        """Method to export the pharmacophore to the pharmagist compatible format.

            Parameters
            ----------
            file_name : str
                Name of file to be written with the pharmagist format of the pharmacophore.
        """
        return to_pharmagist(self, file_name=file_name)

    def to_moe(self, file_name: str) -> None:
        """Method to export the pharmacophore to the MOE compatible format.

            Parameters
            ----------
            file_name: str
                Name of file to be written with the MOE format of the pharmacophore.

        """
        ph4_str = _moe_ph4_string(self._pharmacophoric_points)
        with open(file_name, "w") as f:
            f.writelines(ph4_str)

    def to_rdkit(self) -> Tuple[rdkitPharmacophore.Pharmacophore, List[float]]:
        """ Returns a rdkit pharmacophore with the pharmacophoric_points from the original pharmacophore.
            
            rdkit pharmacophores do not store the pharmacophoric_points radii, so they are returned as well.

            Returns
            -------
            rdkit_pharmacophore : rdkit.Chem.Pharm3D.Pharmacophore
                The rdkit pharmacophore.

            radii : list of float
                List with the radius in angstroms of each pharmacophoric point.
        """
        rdkit_element_name = {  # dictionary to map openpharmacophore feature names to rdkit feature names
            "aromatic ring": "Aromatic",
            "hydrophobicity": "Hydrophobe",
            "hb acceptor": "Acceptor",
            "hb donor": "Donor",
            "positive charge": "PosIonizable",
            "negative charge": "NegIonizable",
        }

        points = []
        radii = []

        for element in self._pharmacophoric_points:
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

    def to_nx_graph(self, d_min: float = 2.0, d_max: float = 20.0, bin_size: float = 1.0) -> nx.Graph:
        """ Obtain a networkx graph representation of the pharmacophore.
        
            The pharmacophore graph is a graph whose nodes are pharmacophoric features and
            its edges are the euclidean distance between those features. The distance is 
            discretized into bins so more molecules can match the pharmacophore.
        
            Parameters
            ----------
            d_min : float
                The minimun distance in angstroms from which two pharmacophoric points are considered different.
            
            d_max : flaot
                The maximum distance in angstroms between pharmacohoric points.
                
            bin_size : float
                The size of the bins that will be used to bin the distances.
        
            Returns
            -------
            pharmacophore_graph : networkx.Graph
                The pharmacophore graph
        """
        pharmacophore_graph = nx.Graph()
        bins = np.arange(d_min, d_max, bin_size)
        # We keep track of feature counts to avoid repeated nodes
        feat_count = {
            "A": 0,
            "D": 0,
            "R": 0,
            "H": 0,
            "P": 0,
            "N": 0,
            "E": 0,
            "I": 0,
        }
        nodes = [""] * len(self)
        for ii, point in enumerate(self):
            name = point.short_name
            feat_count[name] += 1
            nodes[ii] = point.short_name + str(feat_count[name])

        for ii in range(len(self)):
            for jj in range(ii + 1, len(self)):
                distance = distance_between_pharmacophoric_points(self[ii], self[jj])
                binned_distance = discretize(distance, bins)
                pharmacophore_graph.add_edge(nodes[ii],
                                             nodes[jj],
                                             dis=binned_distance)

        assert pharmacophore_graph.number_of_nodes() == len(self.pharmacophoric_points), \
            f"Graph has {pharmacophore_graph.number_of_nodes()} nodes"
        return pharmacophore_graph

    def distance_matrix(self) -> np.ndarray:
        """ Compute the distance matrix of the pharmacophore.
        
            Returns
            -------
            dis_matrix : np.ndarray of shape(n_pharmacophoric_points, n_pharmacophoric_points)
                The distance matrix.
        """
        n_pharmacophoric_points = self.n_pharmacophoric_points
        dis_matrix = np.zeros((n_pharmacophoric_points, n_pharmacophoric_points))

        for ii in range(n_pharmacophoric_points):
            for jj in range(ii, n_pharmacophoric_points):
                if ii == jj:
                    dis_matrix[ii, jj] = 0
                else:
                    distance = distance_between_pharmacophoric_points(
                        self._pharmacophoric_points[ii],
                        self._pharmacophoric_points[jj])
                    dis_matrix[ii, jj] = distance
                    dis_matrix[jj, ii] = distance

        return dis_matrix

    def feature_count(self) -> Dict[str, int]:
        """ Count the number of features ot the same kind in the pharmacophore.
        
            Returns
            -------
            counter : dict
                Dictionary with the count of each feature
        """
        counter = {
            "aromatic ring": 0,
            "hydrophobicity": 0,
            "hb acceptor": 0,
            "hb donor": 0,
            "positive charge": 0,
            "negative charge": 0,
        }

        for element in self._pharmacophoric_points:
            counter[element.feature_name] += 1

        return counter

    def __len__(self) -> int:
        return len(self._pharmacophoric_points)

    def __iter__(self) -> "Pharmacophore":
        self.count = -1
        return self

    def __next__(self) -> PharmacophoricPoint:
        self.count += 1
        if self.count >= self.n_pharmacophoric_points:
            raise StopIteration
        return self._pharmacophoric_points[self.count]

    def __getitem__(self, index: int) -> PharmacophoricPoint:
        return self._pharmacophoric_points[index]

    def __eq__(self, other: "Pharmacophore") -> bool:
        """ Check equality between pharmacophores.

            Assumes that pharmacophoric points are sorted equally in both pharmacophores.
        """
        if isinstance(other, type(self)):
            if self.n_pharmacophoric_points == other.n_pharmacophoric_points:
                for ii in range(self.n_pharmacophoric_points):
                    if not self._pharmacophoric_points[ii] == other._pharmacophoric_points[ii]:
                        return False
                return True
        return False

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n_pharmacophoric_points: {self.n_pharmacophoric_points})"
