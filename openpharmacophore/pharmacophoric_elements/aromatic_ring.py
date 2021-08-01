from .features import AromaticRing
from .shapes import Point, Sphere, SphereAndVector, Shapelet, GaussianKernel

class AromaticRingPoint(AromaticRing, Point):

    def __init__(self, position):

        AromaticRing.__init__(self)
        Point.__init__(self, position)

class AromaticRingSphere(AromaticRing, Sphere):

    def __init__(self, center, radius):

        AromaticRing.__init__(self)
        Sphere.__init__(self, center, radius)

class AromaticRingSphereAndVector(AromaticRing, SphereAndVector):

    def __init__(self, center, radius, direction):

        AromaticRing.__init__(self)
        SphereAndVector.__init__(self, center, radius, direction)

class AromaticRingGaussianKernel(AromaticRing, GaussianKernel):

    def __init__(self, center, sigma):
        AromaticRing.__init__(self)
        GaussianKernel.__init__(self, center, sigma)

class AromaticRingShapelet(AromaticRing, Shapelet):

    def __init__(self):

        AromaticRing.__init__(self)
        Shapelet.__init__(self)

