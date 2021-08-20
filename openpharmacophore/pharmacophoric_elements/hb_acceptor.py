from .features import HBAcceptor
from .shapes import Point, Sphere, SphereAndVector, GaussianKernel

class HBAcceptorPoint(HBAcceptor, Point):

    def __init__(self, position):

        HBAcceptor.__init__(self)
        Point.__init__(self, position)

    def __repr__(self):
        return f"{self.__class__.__name__}(position: {self.position})"

class HBAcceptorSphere(HBAcceptor, Sphere):

    def __init__(self, center, radius):

        HBAcceptor.__init__(self)
        Sphere.__init__(self, center, radius)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; radius: {self.radius})"

class HBAcceptorSphereAndVector(HBAcceptor, SphereAndVector):

    def __init__(self, center, radius, direction):

        HBAcceptor.__init__(self)
        SphereAndVector.__init__(self, center, radius, direction)

    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; radius: {self.radius}; direction: {self.direction})"

class HBAcceptorGaussianKernel(HBAcceptor, GaussianKernel):

    def __init__(self, center, sigma):

        HBAcceptor.__init__(self)
        GaussianKernel.__init__(self, center, sigma)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; sigma: {self.sigma})"

