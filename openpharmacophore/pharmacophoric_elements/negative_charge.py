from .features import NegativeCharge
from .shapes import Point, Sphere, GaussianKernel
import numpy as np
import pyunitwizard as puw

class NegativeChargePoint(NegativeCharge, Point):

    def __init__(self, position):

        NegativeCharge.__init__(self)
        Point.__init__(self, position)
    
    def __eq__(self, other):
        if isinstance(other, type(self)):
            return np.all(self.position == other.position)
        return False
    
    def __repr__(self):
        position = np.around(puw.get_value(self.position, "angstroms"), 2)
        return f"{self.__class__.__name__}(position: {position})"

class NegativeChargeSphere(NegativeCharge, Sphere):

    def __init__(self, center, radius):

        NegativeCharge.__init__(self)
        Sphere.__init__(self, center, radius)
    
    def __eq__(self, other):
        if isinstance(other, type(self)):
            radius_eq = self.radius == other.radius
            center_eq = np.all(self.center == other.center)
            return radius_eq and center_eq
        return False
    

    def __repr__(self):
        # String representation of the class is always in angstroms and rounded to 2 decimals
        center = np.around(puw.get_value(self.center, "angstroms"), 2)
        radius = np.around(puw.get_value(self.radius, "angstroms"), 2)
        return f"{self.__class__.__name__}(center: {center}; radius: {radius})"


class NegativeChargeGaussianKernel(NegativeCharge, GaussianKernel):

    def __init__(self, center, sigma):

        NegativeCharge.__init__(self)
        GaussianKernel.__init__(self, center, sigma)
    
    def __eq__(self, other):
        if isinstance(other, type(self)):
            sigma_eq = self.sigma == other.sigma
            center_eq = np.all(self.center == other.center)
            return sigma_eq and center_eq
        return False

    def __repr__(self):
        center = np.around(puw.get_value(self.center, "angstroms"), 2)
        sigma = np.around(puw.get_value(self.sigma, "angstroms"), 2)
        return f"{self.__class__.__name__}(center: {center}; sigma: {sigma})"

