from openpharmacophore import puw
from openpharmacophore._private_tools.exceptions import *

class Point():

    """ Parent class of pharmacophoric point.

    Common attributes and methods will be included here to be inherited by the specific pharmacophoric
    points classes.

    Parameters
    ----------
    position : Quantity (shape:[3,], dimensionality:[L], value:list, tuple, ndarray)
        Coordinates to set the point position in the three dimensional space.
    feature : str
        Pharmacophoric feature name.

    Attributes
    ----------
    position : Quantity (shape:[3,], dimensionality:[L], value:ndarray) or None
        Coordinates of the point in the three dimensional space.
    feature : str
        Pharmacophoric feature name.

    """

    def __init__(self, position, feature):

        #: The arguments checking should be included with decorators in the future
        if not puw.check(position, dimensionality=['L'], value=[list, tuple, ndarray]):
            raise InputArgumentError('position')
        if not type(feature)==str:
            raise InputArgumentError('feature')

        self.position = puw.standardize(position)
        self.feature = feature

    def add_to_NGLView(self, view, color_palette='openpharmacophore', color=None):
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

        if color is None:
            color = get_color_from_palette_for_point('positive', color_palette)

        radius = 0.05
        center = puw.get_value(self.position, unit='nm').to_list()
        name = self.feature
        view.shape.add_sphere(center, color, radius, name)

        pass

