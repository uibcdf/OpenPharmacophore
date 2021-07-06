from openpharmacophore import puw
from openpharmacophore._private_tools.exceptions import *
from openpharmacophore.pharmacophoric_features.color_palettes import get_color_from_palette_for_feature

class Sphere():

    """ Parent class of pharmacophoric point.

    Common attributes and methods will be included here to be inherited by the specific pharmacophoric
    points classes.

    Parameters
    ----------
    center : Quantity (dimensionality:[L]; value_type:list,tuple,ndarray; shape:[3,])
        Coordinates of the sphere center.
    radius : Quantity (dimensionality:[L]; value:float)
        Radius of the pharmacophoric sphere.
    feature : str
        Pharmacophoric feature name.

    Attributes
    ----------
    center : Quantity (dimensionality:[L]; value:ndarray; shape:[3,]) or None
        Coordinates of the sphere center.
    radius : Quantity (dimensionality:[L]; value:float)
        Radius of the pharmacophoric sphere.
    feature : str
        Pharmacophoric feature name.

    """

    def __init__(self, center, radius):

        #: The arguments checking should be included with decorators in the future
        if not puw.check(center, dimensionality=['L'], value=[list, tuple, ndarray]):
            raise InputArgumentError('center')
        if not puw.check(radius, dimensionality=['L'], value=[float]):
            raise InputArgumentError('radius')

        self.center = puw.standardize(center)
        self.radius = puw.standardize(radius)

    def add_to_NGLView(self, view, color_palette='openpharmacophore', color=None, opacity=0.5):
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

        if color is None:
            color = get_color_from_palette_for_point('positive', color_palette)

        radius = puw.get_value(self.radius, unit='nm')
        center = puw.get_value(self.position, unit='nm').to_list()
        name = self.feature

        view.shape.add_sphere(center, color, radius, name)
        component_index = view.n_components-1
        view.update_representation(component=component_index, repr_index=0, opacity=opacity)

        pass

