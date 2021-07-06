from ._point import Point
from ._sphere import Sphere

class PositiveChargePoint(Point):

    def __init__(self, position):

        super().__init__(position, 'positive charge')

class PositiveChargeSphere(Sphere):

    def __init__(self, center, radius):

        super().__init__(center, radius, 'positive charge')

