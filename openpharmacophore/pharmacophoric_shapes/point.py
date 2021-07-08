import numpy as np
from openpharmacophore import _puw
from openpharmacophore import __documentation_web__
from uibcdf_stdlib.input_arguments import check_input_argument
from uibcdf_stdlib.exceptions import InputArgumentError
from openpharmacophore._private_tools.exceptions import ShapeWithNoColorError
from openpharmacophore.pharmacophoric_features.color_palettes import get_color_from_palette_for_feature

class Point():

    """ Parent class of pharmacophoric point.

    Common attributes and methods will be included here to be inherited by the specific pharmacophoric
    points classes.

    Parameters
    ----------
    position : Quantity (shape:[3,], dimensionality:[L], value:list, tuple, numpy.ndarray)
        Coordinates to set the point position in the three dimensional space.

    Attributes
    ----------
    position : Quantity (shape:[3,], dimensionality:[L], value:numpy.ndarray) or None
        Coordinates of the point in the three dimensional space.

    """

    def __init__(self, position):

        #: The arguments checking should be included with decorators in the future
        #: InputArgumentError shouldn't need arguments
        if not check_input_argument(position, 'quantity', dimensionality=['L'], value_type=[list, tuple, np.ndarray]):
            raise InputArgumentError('position', 'Point', __documentation_web__)

        self.position = _puw.standardize(position)

    def add_to_NGLView(self, view, feature_name=None, color_palette='openpharmacophore', color=None):
        """Adding the point representation to a NGLview view

        Note
        ----
        This method does not return a new view but modifies the input object.

        Parameters
        ----------
        view : NGLView.view object
            NGLview object where the point representations is added.
        color_palette : str or dict, default: 'openpharmacophore'
            Color palette to show the point representation.
        color : str or list
            Color to show the point representation as HEX or RGB code.

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

        radius = 0.05
        center = _puw.get_value(self.position, unit='nm').to_list()

        view.shape.add_sphere(center, color, radius, feature_name)

        pass

