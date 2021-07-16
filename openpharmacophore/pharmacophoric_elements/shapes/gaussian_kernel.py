"""Parent class for pharmacophoric elements with the shape: Gaussian kernel.

This module contains a parent class to be inherited with attributes and methods for pharamacophoric
elements with the 'gaussian kernel' shape.

"""

import numpy as np
from openpharmacophore import _puw
from openpharmacophore import __documentation_web__
from uibcdf_stdlib.input_arguments import check_input_argument
from uibcdf_stdlib.exceptions import InputArgumentError
from uibcdf_stdlib.colors import convert as convert_color_code
from openpharmacophore._private_tools.exceptions import ShapeWithNoColorError
from openpharmacophore.pharmacophoric_elements.features.color_palettes import get_color_from_palette_for_feature

class GaussianKernel():

    """ Parent class for the pharmacophoric shape Gaussian kernel.

    Common attributes and methods will be included here to be inherited by the specific pharmacophoric
    elements with shape Gaussian kernel.

    Parameters
    ----------
    center : Quantity (dimensionality:{'[L]':1}; value_type:list, tuple, numpy.ndarray; shape:(3,))
        Coordinates of the Gaussian kernel center.
    sigma : Quantity (dimensionality:; value:float)
        Standard deviation of the Gaussian kernel of the pharmacophoric sphere.

    Attributes
    ----------
    center : Quantity (dimensionality:{'[L]':1}; value_type:numpy.ndarray; shape:(3,))
        Coordinates of the Gaussian kernel center.
    sigma : Quantity (dimensionality:; value:float)
        Standard deviation of the Gaussian kernel of the pharmacophoric sphere.

    """

    def __init__(self, center, sigma):

        #: The arguments checking should be included with decorators in the future
        #: And the error should probably be raised in the method 'check_input_argument'
        if not check_input_argument(center, 'quantity', dimensionality={'[L]':1}, value_type=[list, tuple, np.ndarray]):
            raise InputArgumentError('center', 'GaussianKernel', __documentation_doc__)
        if not check_input_argument(sigma, 'quantity', dimensionality={'[L]':1}, value_type=float):
            raise InputArgumentError('sigma', 'GaussianKernel', __documentation_doc__)

        self.shape_name = 'gaussian kernel'

        self.center = _puw.standardize(center)
        self.sigma = _puw.standardize(sigma)

    def add_to_NGLView(self, view, feature_name=None, color_palette='openpharmacophore', color=None, opacity=0.5):
        """Adding the Gaussian kernel representation to an NGLview view

        Parameters
        ----------
        view : NGLView.view object
            NGLview object where the point representations is added.
        color_palette : str or dict, default: 'openpharmacophore'
            Color palette to show the Gaussian kernel representation.
        color : str or list
            Color to show the Gaussian kernel representation as HEX or RGB code.

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

        center = _puw.get_value(self.center, to_unit='angstroms').tolist()
        sigma = _puw.get_value(self.sigma, to_unit='angstroms')

        #A Gaussian kernel may be represented as three concentric transparent spheres with radius 0.5 sigma, 1 sigma and 2 sigma

        try:
            n_components = len(view._ngl_component_ids)
        except:
            n_components = 0

        view.shape.add_sphere(center, color, 0.5*sigma)
        view.update_representation(component=n_components, repr_index=0, opacity=opacity)

        view.shape.add_sphere(center, color, sigma)
        view.update_representation(component=n_components+1, repr_index=0, opacity=opacity)

        view.shape.add_sphere(center, color, 2.0*sigma, feature_name)
        view.update_representation(component=n_components+2, repr_index=0, opacity=opacity)

        pass

