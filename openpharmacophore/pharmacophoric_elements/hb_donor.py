from .features import HBDonor
from .shapes import Point, Sphere, SphereAndVector, GaussianKernel

class HBDonorPoint(HBDonor, Point):

    def __init__(self, position):

        HBDonor.__init__(self)
        Point.__init__(self, position)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(position: {self.position})"

class HBDonorSphere(HBDonor, Sphere):

    def __init__(self, center, radius):

        HBDonor.__init__(self)
        Sphere.__init__(self, center, radius)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; radius: {self.radius})"

class HBDonorSphereAndVector(HBDonor, SphereAndVector):

    def __init__(self, center, radius, direction):

        HBDonor.__init__(self)
        SphereAndVector.__init__(self, center, radius, direction)
    
    def __str__(self):
        return f"{self.__class__.__name__}(center: {self.center}; radius: {self.radius}; direction: {self.direction})"

class HBDonorGaussianKernel(HBDonor, GaussianKernel):

    def __init__(self, center, sigma):

        HBDonor.__init__(self)
        GaussianKernel.__init__(self, center, sigma)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; sigma: {self.sigma})"

