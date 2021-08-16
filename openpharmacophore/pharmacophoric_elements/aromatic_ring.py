from .features import AromaticRing
from .shapes import Point, Sphere, SphereAndVector, Shapelet, GaussianKernel

class AromaticRingPoint(AromaticRing, Point):

    def __init__(self, position):

        AromaticRing.__init__(self)
        Point.__init__(self, position)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(position: {self.position})"

class AromaticRingSphere(AromaticRing, Sphere):

    def __init__(self, center, radius):

        AromaticRing.__init__(self)
        Sphere.__init__(self, center, radius)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; radius: {self.radius})"

class AromaticRingSphereAndVector(AromaticRing, SphereAndVector):

    def __init__(self, center, radius, direction):

        AromaticRing.__init__(self)
        SphereAndVector.__init__(self, center, radius, direction)

    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; radius: {self.radius}; direction: {self.direction})"

class AromaticRingGaussianKernel(AromaticRing, GaussianKernel):

    def __init__(self, center, sigma):
        AromaticRing.__init__(self)
        GaussianKernel.__init__(self, center, sigma)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; sigma: {self.sigma})"

class AromaticRingShapelet(AromaticRing, Shapelet):

    def __init__(self):

        AromaticRing.__init__(self)
        Shapelet.__init__(self)
