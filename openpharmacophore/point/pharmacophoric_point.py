# OpenPharmacophore
import openpharmacophore._private_tools.exceptions as exc
from openpharmacophore._private_tools.input_arguments import validate_input_array_like, validate_input_quantity
# Third Party
from matplotlib.colors import to_rgb
import numpy as np
import pyunitwizard as puw


class PharmacophoricPoint:
    """ Class to store pharmacophoric points of any feature type.

        Parameters
        -----------

        feat_type : str
            The feature type of the point. This can be hb donor, hb acceptor, aromatic ring
            hydrophobic, positive charge, negative charge, excluded volume.

        center : Quantity (dimensionality:{'[L]':1}; shape:(3,))
            Coordinates of the sphere center.
        
        radius : Quantity (dimensionality:{'[L]':1})
            Radius of the pharmacophoric sphere.
        
        direction : array-like; shape:(3,), optional
            Direction as a three-dimensional vector.
        
        atom_indices : array-like of int, optional
            The indices of the atoms corresponding to the pharmacophoric point in the molecule from which
            they were extracted.

        Attributes
        ----------
        feature_name : str
            The feature type of the point.

        center : Quantity (dimensionality:{'[L]':1}; shape:(3,))
            Coordinates of the sphere center.

        radius : Quantity (dimensionality:{'[L]':1}; value:float)
            Radius of the pharmacophoric sphere.

        direction : ndarray; shape:(3,)
            Unit vector.
        
        atom_indices : set of int
            A set of the indices of the atoms corresponding to the pharmacophoric point in the molecule from which
            they were extracted.
    
    """
    # TODO : move this dict to another file
    feature_to_char = {
        "hb acceptor": "A",
        "hb donor": "D",
        "aromatic ring": "R",
        "hydrophobicity": "H",
        "positive charge": "P",
        "negative charge": "N",
        "excluded volume": "E",
        "included volume": "I",
    }

    char_to_feature = {
        char: name for name, char in feature_to_char.items()
    }

    palette = {
        'positive charge': '#3498DB',  # Blue
        'negative charge': '#884EA0',  # Purple
        'hb acceptor': '#B03A2E',  # Red
        'hb donor': '#17A589',  # Green
        'included volume': '#707B7C',  # Gray
        'excluded volume': '#283747',  # Black
        'hydrophobicity': '#F5B041',  # Orange
        'aromatic ring': '#F1C40F',  # Yellow
    }

    def __init__(self, feat_type, center, radius, direction=None, atom_indices=None):

        try:
            self._feat_name = PharmacophoricPoint.feature_to_char[feat_type]
        except KeyError:
            raise exc.InvalidFeatureError(feat_type)

        self._validate_input_arguments(center, radius, direction)
        self._center = puw.standardize(center)
        self._radius = puw.standardize(radius)

        if direction is not None:
            self._direction = direction / np.linalg.norm(direction)
        else:
            self._direction = None

        if atom_indices is not None:
            self._atom_indices = set(atom_indices)
        else:
            self._atom_indices = set()

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, new_center):
        validate_input_quantity(new_center, {"[L]": 1}, "center", shape=(3,))
        self._center = new_center

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, new_radius):
        validate_input_quantity(new_radius, {"[L]": 1}, "radius")
        self._radius = new_radius

    @property
    def direction(self):
        return self._direction

    @property
    def atom_indices(self):
        return self._atom_indices

    @atom_indices.setter
    def atom_indices(self, atom_indices):
        self._atom_indices = set(atom_indices)

    @property
    def has_direction(self):
        return self._direction is not None

    @property
    def feature_name(self):
        return self.char_to_feature[self._feat_name]

    @property
    def short_name(self):
        return self._feat_name

    @staticmethod
    def _validate_input_arguments(center, radius, direction):
        """ Validates the input arguments of the PharmacophoricPoint constructor.
        """
        # Validate center
        validate_input_quantity(center, {"[L]": 1}, "center", shape=(3,))
        # Validate radius
        validate_input_quantity(radius, {"[L]": 1}, "radius")
        # Validate direction
        if direction is not None:
            validate_input_array_like(direction, shape=(3,), name="direction")

    def add_to_ngl_view(self, view, color_palette=None, opacity=0.5):
        """ Add the pharmacophoric point to an NGLview.

        Parameters
        ----------
        view : NGLView.widget
            View object where the point representations is added.

        color_palette : dict, default=None
            Color palette to show the point representation.

        opacity : float
            The level of opacity. Must be a number between 0 and 1.

        Notes
        -----
        This method does not return a new view but modifies the input object.

        """
        if color_palette:
            color = color_palette[self.feature_name]
        else:
            color = PharmacophoricPoint.palette[self.feature_name]

        color = to_rgb(color)
        radius = puw.get_value(self._radius, to_unit='angstroms')
        center = puw.get_value(self._center, to_unit='angstroms').tolist()

        n_components = len(view._ngl_component_ids)

        view.shape.add_sphere(center, color, radius, self.feature_name)
        view.update_representation(component=n_components, repr_index=0, opacity=opacity)

        if self.has_direction:
            arrow_radius = 0.2
            end_arrow = puw.get_value(self._center + self._radius * self._direction*2.0,
                                      to_unit='angstroms').tolist()

            view.shape.add_arrow(center, end_arrow, color, arrow_radius)
            view.update_representation(component=n_components + 1, repr_index=0, opacity=0.9)

    @staticmethod
    def get_valid_features():
        """ Get a list of all valid chemical features for a PharmacophoricPoint object"""
        return [
            "hb acceptor",
            "hb donor",
            "aromatic ring",
            "hydrophobicity",
            "positive charge",
            "negative charge",
            "excluded volume",
            "included volume",
        ]

    def is_equal(self, other):
        """ Compare equality of two pharmacophoric points based on atoms indices.
        """
        if self._feat_name == other.short_name and self._atom_indices == other.atom_indices:
            return True
        return False

    def __eq__(self, other):
        """ Compare equality of two pharmacophoric points based on their 3D coordinates,
            directionality and radius.
        """
        if isinstance(other, type(self)):
            if self._feat_name != other._feat_name:
                return False
            if self.has_direction != other.has_direction:
                return False
            radius_eq = np.allclose(self._radius, other.radius, rtol=0, atol=1e-02)
            center_eq = np.allclose(self._center, other.center, rtol=0, atol=1e-04)
            if self.has_direction:
                direction_eq = np.allclose(self._direction, other.direction, rtol=0, atol=1e-04)
                return radius_eq and center_eq and direction_eq
            else:
                return radius_eq and center_eq
        return False

    def __repr__(self):
        center = np.around(puw.get_value(self._center, "angstroms"), 2)
        radius = np.around(puw.get_value(self._radius, "angstroms"), 2)
        x, y, z = center[0], center[1], center[2]
        if self.has_direction:
            direction = np.around(self._direction, 2)
            xd, yd, zd = direction[0], direction[1], direction[2]
            return (f"{self.__class__.__name__}("
                    f"feat_type={self.feature_name}; "
                    f"center=({x}, {y}, {z}); radius={radius}; "
                    f"direction=({xd}, {yd}, {zd}))")
        else:
            return (f"{self.__class__.__name__}("
                    f"feat_type={self.feature_name}; "
                    f"center=({x}, {y}, {z}); radius={radius})")


def distance_between_pharmacophoric_points(p1: PharmacophoricPoint, p2: PharmacophoricPoint):
    """ Compute the distance in angstroms between two pharmacophoric points.
    
        Parameters
        ----------
        p1 : openpharmacophore.pharmacophoric_point
            A pharmacophoric point

        p2 : openpharmacophore.pharmacophoric_point
            A pharmacophoric point
        
        Returns
        -------
        float
            The distance between the pharmacophoric points
    """

    p1_center = puw.get_value(p1.center, "angstroms")
    p2_center = puw.get_value(p2.center, "angstroms")

    return np.linalg.norm(p1_center - p2_center)
