# OpenPharmacophore
from openpharmacophore.pharmacophore.color_palettes import get_color_from_palette_for_feature
from openpharmacophore._private_tools.exceptions import InvalidFeatureType, OpenPharmacophoreTypeError, \
    NegativeRadiusError
from openpharmacophore._private_tools.colors import convert as convert_color_code
from openpharmacophore._private_tools.input_arguments import validate_input_array_like, validate_input_quantity
from openpharmacophore._private_tools.type_vars import ArrayLike, Quantity, Optional, List, Sequence
# Third Party
import numpy as np
import pyunitwizard as puw


class PharmacophoricPoint:
    """ Class to store pharmacophoric points of any feature type. This class can
        store spheres or sphere and vector pharmacophoric points.

        Parameters
        -----------

        feat_type : str
            The feature type of the point. This can be hb donor, hb acceptor, aromatic ring
            hydrophobic, positive charge, negative charge, excluded volume.

        center : Quantity (dimensionality:{'[L]':1}; value_type:list,tuple,numpy.ndarray; shape:(3,))
            Coordinates of the sphere center.
        
        radius : Quantity (dimensionality:{'[L]':1}; value:float)
            Radius of the pharmacophoric sphere.
        
        direction : array_like; shape:(3,)
            Direction as a three-dimensional vector. If the pharmacophoric point doesn't have
            direction is set to None.
        
        atom_indices : list, set or tuple of int
            The indices of the atoms corresponding to the pharmacophoric point in the molecule from which
            they were extracted. A list, set or tuple.

        Attributes
        ----------
        feature_name : str
            The feature type of the point. This can be hb donor, hb acceptor, aromatic ring
            hydrophobic, positive charge, negative charge, excluded volume.

        center : Quantity (dimensionality:{'[L]':1}; value:numpy.ndarray; shape:(3,)) or None
            Coordinates of the sphere center.

        radius : Quantity (dimensionality:{'[L]':1}; value:float)
            Radius of the pharmacophoric sphere.

        direction : ndarray; shape:(3,)
            Unit vector.
        
        atom_indices : set of int
            A set of the indices of the atoms corresponding to the pharmacophoric point in the molecule from which
            they were extracted.
    
    """

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

    def __init__(self, feat_type: str, center: Quantity, radius: Quantity,
                 direction: Optional[ArrayLike] = None, atom_indices: Optional[Sequence] = None) -> None:

        self._validate_input_arguments(feat_type, center, radius, direction, atom_indices)

        self._center = puw.standardize(center)
        self._radius = puw.standardize(radius)
        self._feat_name = PharmacophoricPoint.feature_to_char[feat_type]

        if direction is not None:
            self._direction = direction / np.linalg.norm(direction)
        else:
            self._direction = None

        if atom_indices is not None:
            self._atom_indices = set(atom_indices)
        else:
            self._atom_indices = set()

        self.pharmacophore_index = 0
        self.count = 1

    @property
    def center(self) -> Quantity:
        return self._center

    @property
    def radius(self) -> Quantity:
        return self._radius

    @property
    def direction(self) -> np.ndarray:
        return self._direction

    @property
    def atom_indices(self) -> set:
        return self._atom_indices

    @atom_indices.setter
    def atom_indices(self, atom_indices: Sequence) -> None:
        if not isinstance(atom_indices, (list, set, tuple)):
            raise OpenPharmacophoreTypeError("atom_indices must be a list, set or tuple of int")
        self._atom_indices = atom_indices

    @property
    def has_direction(self) -> bool:
        return self._direction is not None

    @property
    def feature_name(self) -> str:
        return self.char_to_feature[self._feat_name]

    @property
    def short_name(self) -> str:
        return self._feat_name

    @staticmethod
    def _validate_input_arguments(feat_type: str, center: Quantity, radius: Quantity,
                                  direction: ArrayLike, atom_indices: Sequence) -> None:
        """ Validates the input arguments of the PharmacophoricPoint constructor.
        """
        # Validate center
        validate_input_quantity(center, {"[L]": 1}, "center", shape=(3,))
        # Validate radius
        validate_input_quantity(radius, {"[L]": 1}, "radius")
        if puw.get_value(radius, "angstroms") < 0:
            raise NegativeRadiusError("radius must be a positive quantity")
        # Validate direction
        if direction is not None:
            validate_input_array_like(direction, shape=(3,), name="direction")
        # Validate feat type
        if not isinstance(feat_type, str):
            raise InvalidFeatureType("feat_type must be a string")
        # Validate atom_indices
        if atom_indices is not None:
            if not isinstance(atom_indices, (list, set, tuple)):
                raise OpenPharmacophoreTypeError("atom_indices must be a list, set or tuple of int")

        if feat_type not in list(PharmacophoricPoint.feature_to_char.keys()):
            raise InvalidFeatureType(
                f"{feat_type} is not a valid feature type."
                f"Valid feature names are {list(PharmacophoricPoint.feature_to_char.keys())}")

    def add_to_ngl_view(self, view, color_palette='openpharmacophore', opacity=0.5):
        """ Add the pharmacophoric point to an NGLview.

        Parameters
        ----------
        view : NGLView.widget
            View object where the point representations is added.

        color_palette : str or dict, default='openpharmacophore'
            Color palette to show the point representation.

        opacity : float
            The level of opacity. Must be a number between 0 and 1.

        Notes
        -----
        This method does not return a new view but modifies the input object.

        """
        color = get_color_from_palette_for_feature(self.feature_name, color_palette)
        color = convert_color_code(color, to_form='rgb')

        radius = puw.get_value(self._radius, to_unit='angstroms')
        center = puw.get_value(self._center, to_unit='angstroms').tolist()

        try:
            n_components = len(view._ngl_component_ids)
        except:
            n_components = 0

        view.shape.add_sphere(center, color, radius, self.feature_name)
        view.update_representation(component=n_components, repr_index=0, opacity=opacity)

        if self.has_direction:
            arrow_radius = 0.2
            end_arrow = puw.get_value(self._center + self._radius * self._direction, to_unit='angstroms').tolist()

            view.shape.add_arrow(center, end_arrow, color, arrow_radius)
            view.update_representation(component=n_components + 1, repr_index=0, opacity=0.9)

    @staticmethod
    def get_valid_features() -> List[str]:
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

    def is_equal(self, other: "PharmacophoricPoint") -> bool:
        """ Compare equality of two pharmacophoric points based on atoms indices.
        """
        if self._feat_name == other._feat_name and self._atom_indices == other.atom_indices:
            return True
        return False

    def __eq__(self, other: "PharmacophoricPoint") -> bool:
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

    def __repr__(self) -> str:
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


# TODO: we should try to implement dynophores without the use of this class
class UniquePharmacophoricPoint(PharmacophoricPoint):
    """ A class to keep track of unique pharmacophoric points on a dynophore or in a set
        with multiple pharmacophores.

        Inherits from PharmacophoricPoint.

    """

    def __init__(self, point: PharmacophoricPoint, pharmacophore_index: int) -> None:
        super().__init__(point.feature_name,
                         point.center,
                         point.radius,
                         direction=None,
                         atom_indices=point.atom_indices)
        self.count = 1  # To keep the count of each point when working with multiple pharmacophores.
        self.frequency = 0.0
        self.timesteps = [pharmacophore_index]


def distance_between_pharmacophoric_points(p1: PharmacophoricPoint, p2: PharmacophoricPoint) -> float:
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
