from .features import NegativeCharge
from .shapes import Point, Sphere, GaussianKernel

class NegativeChargePoint(NegativeCharge, Point):

    def __init__(self, position):

        NegativeCharge.__init__(self)
        Point.__init__(self, position)

class NegativeChargeSphere(NegativeCharge, Sphere):

    def __init__(self, center, radius):

        NegativeCharge.__init__(self)
        Sphere.__init__(self, center, radius)

class NegativeChargeGaussianKernel(NegativeCharge, GaussianKernel):

    def __init__(self, center, sigma):

        NegativeCharge.__init__(self)
        GaussianKernel.__init__(self, center, sigma)

