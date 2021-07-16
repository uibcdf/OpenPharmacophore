"""Parent class for pharmacophoric elements with the shape: point.

This module contains a parent class to be inherited with attributes and methods for pharamacophoric
elements with the 'point' shape.

"""

import numpy as np
from openpharmacophore import _puw
from openpharmacophore import __documentation_web__
from uibcdf_stdlib.input_arguments import check_input_argument
from uibcdf_stdlib.exceptions import InputArgumentError
from uibcdf_stdlib.colors import convert as convert_color_code
from openpharmacophore._private_tools.exceptions import ShapeWithNoColorError
from openpharmacophore.pharmacophoric_elements.features.color_palettes import get_color_from_palette_for_feature

class Point():

    """Parent class for the pharmacophoric shape point.

    Common attributes and methods will be included here to be inherited by specific pharmacophoric
    elements with shape point.

    Parameters
    ----------
    position : Quantity (shape:(3,), dimensionality:{'[L]':1}, value:list, tuple, numpy.ndarray)
        Coordinates to set the point position in the three dimensional space.

    Attributes
    ----------
    position : Quantity (shape:(3,), dimensionality:{'[L]':1}, value:numpy.ndarray) or None
        Coordinates of the point in the three dimensional space.

    """

    def __init__(self, position):

        #: The arguments checking should be included with decorators in the future
        #: InputArgumentError shouldn't need arguments
        if not check_input_argument(position, 'quantity', dimensionality={'[L]':1}, value_type=[list, tuple, np.ndarray]):
            raise InputArgumentError('position', 'Point', __documentation_web__)

        self.shape_name = 'point'
        self.position = _puw.standardize(position)

    def add_to_NGLView(self, view, feature_name=None, color_palette='openpharmacophore', color=None, opacity=0.5):
        """Adding the point representation to an NGLview view

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

        radius = 0.05
        center = _puw.get_value(self.position, to_unit='angstroms').tolist()

        try:
            n_components = len(view._ngl_components_ids)
        except:
            n_components = 0

        view.shape.add_sphere(center, color, radius, feature_name)
        view.update_representation(component=n_components, repr_index=0, opacity=opacity)

        pass

