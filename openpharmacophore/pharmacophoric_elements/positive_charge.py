"""Pharmacophoric elements for the feature: positive charge.

This module contains all available classes to create pharmacophoric elements for the feature
'positive charge' with different shapes.

Notes
-----
This classes need to be reviewed. Some of them may be removed in the future if they are useless.

"""

from .features import PositiveCharge
from .shapes import Point, Sphere, GaussianKernel

class PositiveChargePoint(PositiveCharge, Point):

    def __init__(self, position):

        PositiveCharge.__init__(self)
        Point.__init__(self, position)

class PositiveChargeSphere(PositiveCharge, Sphere):

    def __init__(self, center, radius):

        PositiveCharge.__init__(self)
        Sphere.__init__(self, center, radius)

class PositiveChargeGaussianKernel(PositiveCharge, GaussianKernel):

    def __init__(self, center, sigma):

        PositiveCharge.__init__(self)
        GaussianKernel.__init__(self, center, sigma)

