from .features import ExcludedVolume
from .shapes import Point, Sphere, GaussianKernel, Shapelet
import numpy as np
import pyunitwizard as puw

class ExcludedVolumePoint(ExcludedVolume, Point):

    def __init__(self, position):

        ExcludedVolume.__init__(self)
        Point.__init__(self, position)
    
    def __eq__(self, other):
        if isinstance(other, type(self)):
            return np.all(self.position == other.position)
        return False

    def __repr__(self):
        position = np.around(puw.get_value(self.position, "angstroms"), 2)
        return f"{self.__class__.__name__}(position: {position})"

class ExcludedVolumeSphere(ExcludedVolume, Sphere):

    def __init__(self, center, radius):

        ExcludedVolume.__init__(self)
        Sphere.__init__(self, center, radius)
    
    def __eq__(self, other):
        if isinstance(other, type(self)):
            radius_eq = np.allclose(self.radius, other.radius, rtol=0, atol=1e-02)
            center_eq = np.allclose(self.center, other.center, rtol=0, atol=1e-04)
            return radius_eq and center_eq
        return False
        
    def __repr__(self):
        # String representation of the class is always in angstroms and rounded to 4 decimals
        center = np.around(puw.get_value(self.center, "angstroms"), 4)
        radius = np.around(puw.get_value(self.radius, "angstroms"), 2)
        x, y, z = center[0], center[1], center[2]
        return f"{self.__class__.__name__}(center: ({x}, {y}, {z}); radius: {radius})"

class ExcludedVolumeGaussianKernel(ExcludedVolume, GaussianKernel):

    def __init__(self, center, sigma):

        ExcludedVolume.__init__(self)
        GaussianKernel.__init__(self, center, sigma)
    
    def __eq__(self, other):
        if isinstance(other, type(self)):
            sigma_eq = self.sigma == other.sigma
            center_eq = np.allclose(self.center, other.center, rtol=1e-04)
            return sigma_eq and center_eq
        return False
    
    def __repr__(self):
        center = np.around(puw.get_value(self.center, "angstroms"), 2)
        sigma = np.around(puw.get_value(self.sigma, "angstroms"), 2)
        return f"{self.__class__.__name__}(center: {center}; sigma: {sigma})"

class ExcludedVolumeShapelet(ExcludedVolume, Shapelet):

    def __init__(self):

        ExcludedVolume.__init__(self)
        Shapelet.__init__(self)

