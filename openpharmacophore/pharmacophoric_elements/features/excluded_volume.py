"""Parent class for pharmacophoric elements with the feature: excluded volume.

This module contains a parent class to be inherited with attributes and methods for pharamacophoric
elements with the 'excluded volume' feature.

"""

class ExcludedVolume():

    """ Parent class of pharmacophoric feature.

    Common attributes and methods to be inherited by the pharmacophoric elements with the 'excluded
    volume' feature.

    Attributes
    ----------
    feature_name : str
        Feature name: 'excluded volume'.

    """

    def __init__(self):

        self.feature_name = 'excluded volume'

