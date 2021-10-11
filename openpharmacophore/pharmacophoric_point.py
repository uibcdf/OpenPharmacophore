from openpharmacophore import __documentation_web__
from uibcdf_stdlib.colors import convert as convert_color_code
from openpharmacophore._private_tools.exceptions import ShapeWithNoColorError
from openpharmacophore.pharmacophoric_elements.features.color_palettes import get_color_from_palette_for_feature
from uibcdf_stdlib.input_arguments import check_input_argument
from uibcdf_stdlib.exceptions import InputArgumentError
import numpy as np
import pyunitwizard as puw

class PharmacophoricPoint():
    """ Class to store pharmacophoric points of any feature type. This class can
        store sheres or sphere and vector points.

        This class will replace the other classes.

        Parameters:
        ----------

        feat_type: str
            The feature type of the point. This can be hb donor, hb acceptor, aromatic ring
            hydrophobic, positive charge, negative charge, excluded volume.

        center : Quantity (dimensionality:{'[L]':1}; value_type:list,tuple,numpy.ndarray; shape:(3,))
            Coordinates of the sphere center.
        
        radius : Quantity (dimensionality:{'[L]':1}; value:float)
            Radius of the pharmacophoric sphere.
        
        direction : list, tuple, ndarray; shape:(3,)
            Vector direction as a three dimensional vector. If the pharmacophoric point doesn't have
            direction is set to None.
        
        atoms_inxs: tuple of int
            The indices of the atoms corresponding to the pharmacophoic point in the molecule from which
            they were extracted. 

        Attributes
        ----------
        feature_name: str
            The feature type of the point. This can be hb donor, hb acceptor, aromatic ring
            hydrophobic, positive charge, negative charge, excluded volume.

        center : Quantity (dimensionality:{'[L]':1}; value:numpy.ndarray; shape:(3,)) or None
            Coordinates of the sphere center.

        radius : Quantity (dimensionality:{'[L]':1}; value:float)
            Radius of the pharmacophoric sphere.

        direction : list, tuple, ndarray; shape:(3,)
            Unit vector.
        
        has_direction: bool
            True if the point has direction.

        element_name: str
            The name of the element contains the feature type and wheter the point is a sphere or a sphere
            and vector.
        
         atoms_inxs: tuple of int
            The indices of the atoms corresponding to the pharmacophoic point in the molecule from which
            they were extracted. 

    
    """
    def __init__(self, feat_type, center, radius, direction=None, atoms_inxs=None):

        #: InputArgumentError shouldn't need arguments
        if not check_input_argument(center, 'quantity', dimensionality={'[L]':1}, value_type=[list, tuple, np.ndarray]):
            raise InputArgumentError('center', 'SphereAndVector', __documentation_web__)
        if not check_input_argument(radius, 'quantity', dimensionality={'[L]':1}, value_type=[float, int]):
            raise InputArgumentError('radius', 'SphereAndVector', __documentation_web__)
        if not isinstance(feat_type, str):
            raise ValueError("feat_type must be a string")
        
        valid_feat_types = ["aromatic ring",
            "hydrophobicity",
            "hb acceptor",
            "hb donor",
            "positive charge",
            "negative charge",
            "excluded sphere",
            "included sphere"]

        if feat_type not in valid_feat_types:
            raise ValueError(f"{feat_type} is not a valid feature type. Valid feature names are {valid_feat_types}")

        self.center = puw.standardize(center)
        self.radius = puw.standardize(radius)
        self.feature_name = feat_type
        self.element_name = "".join([n.capitalize() for n in self.feature_name.split()])
        if direction is not None:
            if not check_input_argument(direction, [tuple, list, np.ndarray], shape=(3,)):
                raise InputArgumentError('direction', 'SphereAndVector', __documentation_web__)
            self.direction = direction/np.linalg.norm(direction)
            self.has_direction = True
            self.element_name += "SphereAndVector"
        else:
            self.direction = None
            self.has_direction = False
            self.element_name += "Sphere"
        
        if atoms_inxs:
            self.atoms_inxs = atoms_inxs
        else:
            self.atoms_inxs = None
        
    
    def add_to_NGLView(self, view, feature_name=None, color_palette='openpharmacophore', color=None, opacity=0.5):
        """Adding the element representation to an NGLview view

        Parameters
        ----------
        view : NGLView.widget object
            NGLview object where the point representations is added.
        color_palette : str or dict, default: 'openpharmacophore'
            Color palette to show the point representation.
        color : str or list
            Color to show the point representation as HEX or RGB code.

        Note
        ----
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
                raise ShapeWithNoColorError(__documentation_web__)

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

    def get_element_name(self):
        """ Get full element name

            Returns
            -------
            str
                The name of the element
        """
        return self.element_name

    def get_direction(self):
        """ Get the pharmacophoric point direction.

            Returns
            -------
            numpy.ndarray; shape(3,)
                The direction vector.
        """
        if not self.has_direction:
            raise ValueError("This pharmacophoric point has no direction vector")
        return self.direction

    def get_center(self, unit="angstroms"):
        """ Get the pharmacophoric point centroid.

            Returns
            -------
            numpy.ndarray; shape(3,)
                Array with the centroid x, y and z coordinates.
        """
        return puw.get_value(self.center, to_unit=unit)

    def get_radius(self, unit="angstroms"):
        """ Get the pharmacophoric point radius.

            Returns
            -------
            float
                The radius of the pharmacophoric point.
        """
        return puw.get_value(self.radius, to_unit=unit)
    
    def __eq__(self, other):
        if isinstance(other, type(self)):
            if self.feature_name != other.feature_name:
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

    def __repr__(self):
        center = np.around(puw.get_value(self.center, "angstroms"), 4)
        radius = np.around(puw.get_value(self.radius, "angstroms"), 2)
        x, y, z = center[0], center[1], center[2]
        if self.has_direction:
            direction = np.around(self.direction, 4)
            xd, yd, zd = direction[0], direction[1], direction[2]
            return f"{self.element_name}(center: ({x}, {y}, {z}); radius: {radius}; direction: ({xd}, {yd}, {zd}))"
        else:
            return f"{self.element_name}(center: ({x}, {y}, {z}); radius: {radius})"