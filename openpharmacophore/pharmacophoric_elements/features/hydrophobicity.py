"""Parent class for pharmacophoric elements with the feature: positive charge.

This module contains a parent class to be inherited with attributes and methods for pharamacophoric
elements with the 'positive charge' feature.

"""

class Hydrophobicity():

    """ Parent class of pharmacophoric feature.

    Common attributes and methods to be inherited by the pharmacophoric elements with the
    'hydrophobicity' feature.

    Attributes
    ----------
    feature_name : str
        Feature name: 'hydrophobicity'.

    """

    def __init__(self):

        self.feature_name = 'hydrophobicity'

