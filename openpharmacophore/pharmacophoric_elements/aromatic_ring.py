from .features import AromaticRing
from .shapes import Point, Sphere, SphereAndVector, Shapelet, GaussianKernel
import numpy as np
import pyunitwizard as puw


class AromaticRingPoint(AromaticRing, Point):

    def __init__(self, position):

        AromaticRing.__init__(self)
        Point.__init__(self, position)
    
    def __repr__(self):
        position = np.around(puw.get_value(self.position, "angstroms"), 2)
        return f"{self.__class__.__name__}(position: {position})"

class AromaticRingSphere(AromaticRing, Sphere):

    def __init__(self, center, radius):

        AromaticRing.__init__(self)
        Sphere.__init__(self, center, radius)
    
    def __repr__(self):
        # String representation of the class is always in angstroms and rounded to 2 decimals
        center = np.around(puw.get_value(self.center, "angstroms"), 2)
        radius = np.around(puw.get_value(self.radius, "angstroms"), 2)
        return f"{self.__class__.__name__}(center: {center}; radius: {radius})"

class AromaticRingSphereAndVector(AromaticRing, SphereAndVector):

    def __init__(self, center, radius, direction):

        AromaticRing.__init__(self)
        SphereAndVector.__init__(self, center, radius, direction)

    def __repr__(self):
        center = np.around(puw.get_value(self.center, "angstroms"), 2)
        radius = np.around(puw.get_value(self.radius, "angstroms"), 2)
        direction = np.around(self.direction, 2)
        return f"{self.__class__.__name__}(center: {center}; radius: {radius}; direction: {direction})"

class AromaticRingGaussianKernel(AromaticRing, GaussianKernel):

    def __init__(self, center, sigma):
        AromaticRing.__init__(self)
        GaussianKernel.__init__(self, center, sigma)
    
    def __repr__(self):
        center = np.around(puw.get_value(self.center, "angstroms"), 2)
        sigma = np.around(puw.get_value(self.sigma, "angstroms"), 2)
        return f"{self.__class__.__name__}(center: {center}; sigma: {sigma})"

class AromaticRingShapelet(AromaticRing, Shapelet):

    def __init__(self):

        AromaticRing.__init__(self)
        Shapelet.__init__(self)
