from .features import Hydrophobicity
from .shapes import Point, Sphere, GaussianKernel, Shapelet

class HydrophobicPoint(Hydrophobicity, Point):

    def __init__(self, position):

        Hydrophobicity.__init__(self)
        Point.__init__(self, position)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(position: {self.position})"

class HydrophobicSphere(Hydrophobicity, Sphere):

    def __init__(self, center, radius):

        Hydrophobicity.__init__(self)
        Sphere.__init__(self, center, radius)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; radius: {self.radius})"

class HydrophobicGaussianKernel(Hydrophobicity, GaussianKernel):

    def __init__(self, center, sigma):

        Hydrophobicity.__init__(self)
        GaussianKernel.__init__(self, center, sigma)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; sigma: {self.sigma})"

class HydrophobicShapelet(Hydrophobicity, Shapelet):

    def __init__(self):

        Hydrophobicity.__init__(self)
        Shapelet.__init__(self)

