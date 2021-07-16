"""Parent class for pharmacophoric elements with the feature: hb donor.

This module contains a parent class to be inherited with attributes and methods for pharamacophoric
elements with the 'hb donor' feature.

"""

class HBDonor():

    """ Parent class of pharmacophoric feature.

    Common attributes and methods to be inherited by the pharmacophoric elements with the 'hb donor' feature.

    Attributes
    ----------
    feature_name : str
        Feature name: 'hb donor'.

    """

    def __init__(self):

        self.feature_name = 'hb donor'

