from .features import Hydrophobicity
from .shapes import Point, Sphere, GaussianKernel, Shapelet

class HydrophobicPoint(Hydrophobicity, Point):

    def __init__(self, position):

        Hydrophobicity.__init__(self)
        Point.__init__(self, position)

class HydrophobicSphere(Hydrophobicity, Sphere):

    def __init__(self, center, radius):

        Hydrophobicity.__init__(self)
        Sphere.__init__(self, center, radius)

class HydrophobicGaussianKernel(Hydrophobicity, GaussianKernel):

    def __init__(self, center, sigma):

        Hydrophobicity.__init__(self)
        GaussianKernel.__init__(self, center, sigma)

class HydrophobicShapelet(Hydrophobicity, Shapelet):

    def __init__(self):

        Hydrophobicity.__init__(self)
        Shapelet.__init__(self)

