from .features import ExcludedVolume
from .shapes import Point, Sphere, GaussianKernel, Shapelet

class ExcludedVolumePoint(ExcludedVolume, Point):

    def __init__(self, position):

        ExcludedVolume.__init__(self)
        Point.__init__(self, position)

class ExcludedVolumeSphere(ExcludedVolume, Sphere):

    def __init__(self, center, radius):

        ExcludedVolume.__init__(self)
        Sphere.__init__(self, center, radius)

class ExcludedVolumeGaussianKernel(ExcludedVolume, GaussianKernel):

    def __init__(self, center, sigma):

        ExcludedVolume.__init__(self)
        GaussianKernel.__init__(self, center, sigma)

class ExcludedVolumeShapelet(ExcludedVolume, Shapelet):

    def __init__(self):

        ExcludedVolume.__init__(self)
        Shapelet.__init__(self)

