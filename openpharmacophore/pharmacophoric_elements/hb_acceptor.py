from .features import HBAcceptor
from .shapes import Point, Sphere, SphereAndVector, GaussianKernel

class HBAcceptorPoint(HBAcceptor, Point):

    def __init__(self, position):

        HBAcceptor.__init__(self)
        Point.__init__(self, position)

class HBAcceptorSphere(HBAcceptor, Sphere):

    def __init__(self, center, radius):

        HBAcceptor.__init__(self)
        Sphere.__init__(self, center, radius)

class HBAcceptorSphereAndVector(HBAcceptor, SphereAndVector):

    def __init__(self, center, radius, direction):

        HBAcceptor.__init__(self)
        SphereAndVector.__init__(self, center, radius, direction)

class HBAcceptorGaussianKernel(HBAcceptor, GaussianKernel):

    def __init__(self, center, sigma):

        HBAcceptor.__init__(self)
        GaussianKernel.__init__(self, center, sigma)

