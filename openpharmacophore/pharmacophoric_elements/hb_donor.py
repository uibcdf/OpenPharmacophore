from .features import HBDonor
from .shapes import Point, Sphere, SphereAndVector, GaussianKernel

class HBDonorPoint(HBDonor, Point):

    def __init__(self, position):

        HBDonor.__init__(self)
        Point.__init__(self, position)

class HBDonorSphere(HBDonor, Sphere):

    def __init__(self, center, radius):

        HBDonor.__init__(self)
        Sphere.__init__(self, center, radius)

class HBDonorSphereAndVector(HBDonor, SphereAndVector):

    def __init__(self, center, radius, direction):

        HBDonor.__init__(self)
        SphereAndVector.__init__(self, center, radius, direction)

class HBDonorGaussianKernel(HBDonor, GaussianKernel):

    def __init__(self, center, sigma):

        HBDonor.__init__(self)
        GaussianKernel.__init__(self, center, sigma)

