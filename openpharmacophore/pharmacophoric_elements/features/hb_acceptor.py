"""Parent class for pharmacophoric elements with the feature: hb acceptor.

This module contains a parent class to be inherited with attributes and methods for pharamacophoric
elements with the 'hb acceptor' feature.

"""

class HBAcceptor():

    """ Parent class of pharmacophoric feature.

    Common attributes and methods to be inherited by the pharmacophoric elements with the 'hb
    acceptor' feature.

    Attributes
    ----------
    feature_name : str
        Feature name: 'hb acceptor'.

    """

    def __init__(self):

        self.feature_name = 'hb acceptor'

