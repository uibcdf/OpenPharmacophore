# OpenPharmacophore
from openpharmacophore import __documentation_web__
from openpharmacophore.pharmacophore.color_palettes import get_color_from_palette_for_feature
from openpharmacophore._private_tools.exceptions import InvalidFeatureType, OpenPharmacophoreTypeError, PointWithNoColorError, NegativeRadiusError
from openpharmacophore._private_tools.colors import convert as convert_color_code
from openpharmacophore._private_tools.input_arguments import validate_input_array_like, validate_input_quantity
# Third Party
import numpy as np
import pyunitwizard as puw
from pint import Quantity
# Standard library
from typing import List, Optional, TypeVar, Sequence

feature_to_char = {
            "hb acceptor":     "A",
            "hb donor":        "D",
            "aromatic ring":   "R",
            "hydrophobicity":  "H",
            "positive charge": "P",
            "negative charge": "N",
            "excluded volume": "E",
            "included volume": "I",
        }

ArrayLike = TypeVar("ArrayLike", list, np.ndarray)

class PharmacophoricPoint():
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
            Vector direction as a three dimensional vector. If the pharmacophoric point doesn't have
            direction is set to None. It must be of length or shape 3.
        
        atom_indices : list, set or tuple of int
            The indices of the atoms corresponding to the pharmacophoic point in the molecule from which
            they were extracted. A list, set or tuple can be passed. 

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
        
        has_direction : bool
            Whether the pharmacophoric point has direction.
        
        atom_indices : set of int
            A set of the indices of the atoms corresponding to the pharmacophoic point in the molecule from which
            they were extracted.  

    
    """
    def __init__(self, feat_type: str, center: Quantity, radius: Quantity, 
                direction: Optional[ArrayLike] = None, atom_indices: Optional[Sequence] = None) -> None:
        
        # Validate center
        validate_input_quantity(center, {"[L]" : 1}, "center", shape=(3,))
        # Validate radius
        validate_input_quantity(radius, {"[L]" : 1}, "radius")
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
            
        if feat_type not in list(feature_to_char.keys()):
            raise InvalidFeatureType(f"{feat_type} is not a valid feature type. Valid feature names are {list(feature_to_char.keys())}")

        self.center = puw.standardize(center)
        self.radius = puw.standardize(radius)
        self.feature_name = feat_type
        self.short_name = feature_to_char[feat_type]

        if direction is not None:
            self.direction = direction/np.linalg.norm(direction)
            self.has_direction = True
        else:
            self.direction = None
            self.has_direction = False
        
        if atom_indices is not None:
            self.atom_indices = set(atom_indices)
        else:
            self.atom_indices = None
        
        self.pharmacophore_index = 0
        self.count = 1
        
        
    def add_to_NGLView(self, view, feature_name=None, color_palette='openpharmacophore', color=None, opacity=0.5):
        """ Adding the pharmacophoric point to an NGLview view.

        Parameters
        ----------
        view : NGLView.widget
            NGLview object where the point representations is added.
        color_palette : str or dict, default='openpharmacophore'
            Color palette to show the point representation.
        color : str or list
            Color to show the point representation as HEX or RGB code.
        opacity : float
            The level of opacity. Must be a number between 0 and 1.

        Notes
        -----
        This method does not return a new view but modifies the input object.

        """

        if feature_name is None:
            try:
                feature_name = self.feature_name
            except:
                pass

        if color is None:
            if feature_name is not None:
                color = get_color_from_palette_for_feature(feature_name, color_palette)
            else:
                raise PointWithNoColorError(__documentation_web__)

        color = convert_color_code(color, to_form='rgb')

        radius = puw.get_value(self.radius, to_unit='angstroms')
        center = puw.get_value(self.center, to_unit='angstroms').tolist()

        try:
            n_components = len(view._ngl_component_ids)
        except:
            n_components = 0

        view.shape.add_sphere(center, color, radius, feature_name)
        view.update_representation(component=n_components, repr_index=0, opacity=opacity)

        if self.has_direction:
            arrow_radius = 0.2
            end_arrow = puw.get_value(self.center+self.radius*self.direction, to_unit='angstroms').tolist()

            view.shape.add_arrow(center, end_arrow, color, arrow_radius)
            view.update_representation(component=n_components+1, repr_index=0, opacity=0.9)
    
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
        if self.short_name == other.short_name and self.atom_indices == other.atom_indices:
                return True
        return False

    def __eq__(self, other: "PharmacophoricPoint") -> bool:
        """ Compare equality of two pharmacophoric points based on their 3D coordinates,
            directionality and radius.
        """
        if isinstance(other, type(self)):
            if self.short_name != other.short_name:
                return False
            if self.has_direction != other.has_direction:
                return False
            radius_eq = np.allclose(self.radius, other.radius, rtol=0, atol=1e-02)
            center_eq = np.allclose(self.center, other.center, rtol=0, atol=1e-04)
            if self.has_direction:
                direction_eq = np.allclose(self.direction, other.direction, rtol=0, atol=1e-04)
                return radius_eq and center_eq and direction_eq
            else:
                return radius_eq and center_eq
        return False

    def __repr__(self) -> str:
        center = np.around(puw.get_value(self.center, "angstroms"), 2)
        radius = np.around(puw.get_value(self.radius, "angstroms"), 2)
        x, y, z = center[0], center[1], center[2]
        if self.has_direction:
            direction = np.around(self.direction, 4)
            xd, yd, zd = direction[0], direction[1], direction[2]
            return (f"{self.__class__.__name__}("
                    f"feat_type={self.feature_name}; " 
                    f"center=({x}, {y}, {z}); radius: {radius}; " 
                    f"direction=({xd}, {yd}, {zd}))")
        else:
            return (f"{self.__class__.__name__}("
                    f"feat_type={self.feature_name}; " 
                    f"center=({x}, {y}, {z}); radius={radius})")


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
        self.count = 1 # To keep the count of each point when working with multiple pharmacophores.
        self.frequency = 0.0
        self.timesteps = [pharmacophore_index] 


def distance_bewteen_pharmacophoric_points(p1: PharmacophoricPoint, p2: PharmacophoricPoint) -> float:
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