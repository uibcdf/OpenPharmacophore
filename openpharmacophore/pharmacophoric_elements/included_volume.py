from .features import IncludedVolume
from .shapes import Point, Sphere, GaussianKernel, Shapelet

class IncludedVolumePoint(IncludedVolume, Point):

    def __init__(self, position):

        IncludedVolume.__init__(self)
        Point.__init__(self, position)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(position: {self.position})"

class IncludedVolumeSphere(IncludedVolume, Sphere):

    def __init__(self, center, radius):

        IncludedVolume.__init__(self)
        Sphere.__init__(self, center, radius)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; radius: {self.radius})"

class IncludedVolumeGaussianKernel(IncludedVolume, GaussianKernel):

    def __init__(self, center, sigma):

        IncludedVolume.__init__(self)
        GaussianKernel.__init__(self, center, sigma)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}; sigma: {self.sigma})"

class IncludedVolumeShapelet(IncludedVolume, Shapelet):

    def __init__(self):

        IncludedVolume.__init__(self)
        Shapelet.__init__(self)

