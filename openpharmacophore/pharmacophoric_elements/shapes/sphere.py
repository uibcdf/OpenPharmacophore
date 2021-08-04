"""Parent class for pharmacophoric elements with the shape: sphere.

This module contains a parent class to be inherited with attributes and methods for pharamacophoric
elements with the 'sphere' shape.

"""


import numpy as np
from openpharmacophore import _puw
from openpharmacophore import __documentation_web__
from uibcdf_stdlib.input_arguments import check_input_argument
from uibcdf_stdlib.exceptions import InputArgumentError
from uibcdf_stdlib.colors import convert as convert_color_code
from openpharmacophore._private_tools.exceptions import ShapeWithNoColorError
from openpharmacophore.pharmacophoric_elements.features.color_palettes import get_color_from_palette_for_feature

class Sphere():

    """ Parent class for the pharmacophoric shape sphere.

    Common attributes and methods will be included here to be inherited by specific pharmacophoric
    elements with shape sphere.

    Parameters
    ----------
    center : Quantity (dimensionality:{'[L]':1}; value_type:list,tuple,numpy.ndarray; shape:(3,))
        Coordinates of the sphere center.
    radius : Quantity (dimensionality:{'[L]':1}; value:float)
        Radius of the pharmacophoric sphere.

    Attributes
    ----------
    center : Quantity (dimensionality:{'[L]':1}; value:ndarray; shape:(3,)) or None
        Coordinates of the sphere center.
    radius : Quantity (dimensionality:{'[L]':1}; value:float)
        Radius of the pharmacophoric sphere.

    """

    def __init__(self, center, radius):

        #: The arguments checking should be included with decorators in the future
        #: InputArgumentError shouldn't need arguments
        if not check_input_argument(center, 'quantity', dimensionality={'[L]':1}, value_type=[list, tuple, np.ndarray]):
            raise InputArgumentError('center', 'Sphere', __documentation_web__)
        if not check_input_argument(radius, 'quantity', dimensionality={'[L]':1}, value_type=[float, int]):
            raise InputArgumentError('radius', 'Sphere', __documentation_web__)

        self.shape_name = 'sphere'

        self.center = _puw.standardize(center)
        self.radius = _puw.standardize(radius)

    def add_to_NGLView(self, view, feature_name=None, color_palette='openpharmacophore', color=None, opacity=0.5):
        """Adding the sphere representation to an NGLview view

        Parameters
        ----------
        view : NGLView.view object
            NGLview object where the point representations is added.
        color_palette : str or dict, default: 'openpharmacophore'
            Color palette to show the point representation.
        color : str or list
            Color to show the point representation as HEX or RGB code.

        Note
        ----
        This method does not return a new view but modifies in place the input one.

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

        center = _puw.get_value(self.center, to_unit='angstroms').tolist()
        radius = _puw.get_value(self.radius, to_unit='angstroms')

        try:
            n_components = len(view._ngl_component_ids)
        except:
            n_components = 0

        view.shape.add_sphere(center, color, radius, feature_name)
        view.update_representation(component=n_components, repr_index=0, opacity=opacity)

        pass

