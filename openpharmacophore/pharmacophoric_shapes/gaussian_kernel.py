import numpy as np
from openpharmacophore import _puw
from openpharmacophore import __documentation_web__
from uibcdf_stdlib.input_arguments import check_input_argument
from uibcdf_stdlib.exceptions import InputArgumentError
from openpharmacophore._private_tools.exceptions import ShapeWithNoColorError
from openpharmacophore.pharmacophoric_features.color_palettes import get_color_from_palette_for_feature

class GaussianKernel():

    """ Parent class of pharmacophoric Gaussian kernel.

    Common attributes and methods will be included here to be inherited by the specific pharmacophoric
    Gaussian kernels classes.

    Parameters
    ----------
    center : Quantity (dimensionality:[L]; value_type:list, tuple, numpy.ndarray; shape:[3,])
        Coordinates of the Gaussian kernel center.
    sigma : Quantity (dimensionality:[L]; value:float)
        Standard deviation of the Gaussian kernel of the pharmacophoric sphere.

    Attributes
    ----------
    center : Quantity (dimensionality:[L]; value_type:numpy.ndarray; shape:[3,])
        Coordinates of the Gaussian kernel center.
    sigma : Quantity (dimensionality:[L]; value:float)
        Standard deviation of the Gaussian kernel of the pharmacophoric sphere.

    """

    def __init__(self, center, sigma):

        #: The arguments checking should be included with decorators in the future
        #: And the error should probably be raised in the method 'check_input_argument'
        if not check_input_argument(center, 'quantity', dimensionality=['L'], value_type=[list, tuple, np.ndarray]):
            raise InputArgumentError('center', 'GaussianKernel', __documentation_doc__)
        if not check_input_argument(sigma, 'quantity', dimensionality=['L'], value_type=float):
            raise InputArgumentError('sigma', 'GaussianKernel', __documentation_doc__)

        self.center = _puw.standardize(center)
        self.sigma = _puw.standardize(sigma)

    def add_to_NGLView(self, view, feature_name=None, color_palette='openpharmacophore', color=None, opacity=0.5):
        """Adding the Gaussian kernel representation to a NGLview view

        Note
        ----
        This method does not return a new view but modifies the input object.

        Parameters
        ----------
        view : NGLView.view object
            NGLview object where the point representations is added.
        color_palette : str or dict, default: 'openpharmacophore'
            Color palette to show the Gaussian kernel representation.
        color : str or list
            Color to show the Gaussian kernel representation as HEX or RGB code.

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

        opacity_sphere = opacity

        center = _puw.get_value(self.center, unit='nm').to_list()
        sigma = _puw.get_value(self.radius, unit='nm')

        #A Gaussian kernel may be represented as three concentric transparent spheres with radius 0.5 sigma, 1 sigma and 2 sigma

        view.shape.add_sphere(center, color, 0.5*sigma)
        component_index = view.n_components-1
        view.update_representation(component=component_index, repr_index=0, opacity=opacity_sphere)

        view.shape.add_sphere(center, color, sigma)
        component_index = view.n_components-1
        view.update_representation(component=component_index, repr_index=0, opacity=opacity_sphere)

        view.shape.add_sphere(center, color, 2.0*sigma, feature_name)
        component_index = view.n_components-1
        view.update_representation(component=component_index, repr_index=0, opacity=opacity_sphere)

        pass

