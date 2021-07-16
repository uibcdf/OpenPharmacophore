"""Parent class for pharmacophoric elements with the feature: aromatic ring.

This module contains a parent class to be inherited with attributes and methods for pharamacophoric
elements with the 'aromatic ring' feature.

"""

class AromaticRing():

    """ Parent class of pharmacophoric feature.

    Common attributes and methods to be inherited by the pharmacophoric elements with the 'aromatic
    ring' feature.

    Attributes
    ----------
    feature_name : str
        Feature name: 'aromatic ring'.

    """

    def __init__(self):

        self.feature_name = 'aromatic ring'

