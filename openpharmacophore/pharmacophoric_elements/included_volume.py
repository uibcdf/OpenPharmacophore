from openpharmacophore.pharmacophoric_features import IncludedVolume
from openpharmacophore.pharmacophoric_shapes import Point, Sphere, GaussianKernel, Shapelet

class IncludedVolumePoint(IncludedVolume, Point):

    def __init__(self, position):

        IncludedVolume.__init__(self)
        Point.__init__(self, position)

class IncludedVolumeSphere(IncludedVolume, Sphere):

    def __init__(self, center, radius):

        IncludedVolume.__init__(self)
        Sphere.__init__(self, center, radius)

class IncludedVolumeGaussianKernel(IncludedVolume, GaussianKernel):

    def __init__(self, center, sigma):

        IncludedVolume.__init__(self)
        GaussianKernel.__init__(self, center, sigma)

class IncludedVolumeShapelet(IncludedVolume, Shapelet):

    def __init__(self):

        IncludedVolume.__init__(self)
        Shapelet.__init__(self)

