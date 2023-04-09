# OpenPharmacophore
from openpharmacophore.pharmacophore.input_arguments import validate_input_array_like, validate_input_quantity
from openpharmacophore.pharmacophore.exceptions import InvalidFeatureError
from openpharmacophore.constants import PALETTE, FEAT_TO_CHAR, CHAR_TO_FEAT
# Third Party
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
    def __init__(self, feat_type, center, radius, direction=None, atom_indices=None):

        try:
            self._feat_name = FEAT_TO_CHAR[feat_type]
        except KeyError:
            raise InvalidFeatureError(feat_type)

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
        return CHAR_TO_FEAT[self._feat_name]

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

    @staticmethod
    def get_valid_features():
        """ Get a list of all valid chemical features for a PharmacophoricPoint object"""
        return list(PALETTE.keys())

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
