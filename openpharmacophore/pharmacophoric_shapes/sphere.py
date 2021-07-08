import numpy as np
from openpharmacophore import _puw
from openpharmacophore import __documentation_web__
from uibcdf_stdlib.input_arguments import check_input_argument
from uibcdf_stdlib.exceptions import InputArgumentError
from openpharmacophore._private_tools.exceptions import ShapeWithNoColorError
from openpharmacophore.pharmacophoric_features.color_palettes import get_color_from_palette_for_feature

class Sphere():

    """ Parent class of pharmacophoric point.

    Common attributes and methods will be included here to be inherited by the specific pharmacophoric
    points classes.

    Parameters
    ----------
    center : Quantity (dimensionality:[L]; value_type:list,tuple,numpy.ndarray; shape:[3,])
        Coordinates of the sphere center.
    radius : Quantity (dimensionality:[L]; value:float)
        Radius of the pharmacophoric sphere.

    Attributes
    ----------
    center : Quantity (dimensionality:[L]; value:ndarray; shape:[3,]) or None
        Coordinates of the sphere center.
    radius : Quantity (dimensionality:[L]; value:float)
        Radius of the pharmacophoric sphere.

    """

    def __init__(self, center, radius):

        #: The arguments checking should be included with decorators in the future
        #: InputArgumentError shouldn't need arguments
        if not check_input_argument(center, 'quantity', dimensionality=['L'], value_type=[list, tuple, np.ndarray]):
            raise InputArgumentError('center', 'Sphere', __documentation_web___)
        if not check_input_argument(radius, 'quantity', dimensionality=['L'], value_type=float):
            raise InputArgumentError('radius', 'Sphere', __documentation_web__)

        self.center = _puw.standardize(center)
        self.radius = _puw.standardize(radius)

    def add_to_NGLView(self, view, feature_name=None, color_palette='openpharmacophore', color=None, opacity=0.5):
        """Adding the sphere representation to a NGLview view

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

        radius = _puw.get_value(self.radius, unit='nm')
        center = _puw.get_value(self.position, unit='nm').to_list()

        view.shape.add_sphere(center, color, radius, feature_name)
        component_index = view.n_components-1
        view.update_representation(component=component_index, repr_index=0, opacity=opacity)

        pass

