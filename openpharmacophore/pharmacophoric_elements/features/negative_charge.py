
"""Parent class for pharmacophoric elements with the feature: negative charge.

This module contains a parent class to be inherited with attributes and methods for pharamacophoric
elements with the 'negative charge' feature.

"""

class NegativeCharge():

    """ Parent class of pharmacophoric feature.

    Common attributes and methods to be inherited by the pharmacophoric elements with the 'negative
    charge' feature.

    Attributes
    ----------
    feature_name : str
        Feature name: 'negative charge'.

    """

    def __init__(self):

        self.feature_name = 'negative charge'

