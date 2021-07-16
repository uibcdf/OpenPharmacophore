
"""Parent class for pharmacophoric elements with the feature: included volume.

This module contains a parent class to be inherited with attributes and methods for pharmacophoric
elements with the 'included volume' feature.

"""

class IncludedVolume():

    """ Parent class of pharmacophoric feature.

    Common attributes and methods to be inherited by the pharmacophoric elements with the 'included
    volume' feature.

    Attributes
    ----------
    feature_name : str
        Feature name: 'included volume'.

    """

    def __init__(self):

        self.feature_name = 'included volume'

