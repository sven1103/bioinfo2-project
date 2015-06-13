"""
Represents a straight line in 3D space.
"""
import numpy as np


class Line(object):
    """ Represent a straight line in 3D space."""

    def __init__(self, strut, direction):
        """
        Creates straight line in 3D space.

        :param strut: Coordinates of strut vector (x,y,z)
        :type strut: Numeric array (Float0, size 3)
        :param direction: Coordinates of vector with direction of line (x,y,z)
        :type direction: Numeric array (Float0, size 3)
        :return: New instance of Line object.
        """

        self.strut = np.array(strut)
        self.direction = np.array(direction)

    def __str__(self):
        return str(self.strut) + " + " + "a" + str(self.direction)
