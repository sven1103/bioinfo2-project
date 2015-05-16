"""
Implements operations on planes
"""
__author__ = 'lukas'
import numpy as np

class Line(object):
    """ Represent a straight line in 3D space."""

    def __init__(self, strut, direction):

        self.strut = strut
        self.direction = direction

    def __str__(self):
        return str(self.strut) + " + " + "a" + str(self.direction)


class Plane(object):
    """
    Represents plane in 3D space
    """

    def __init__(self, strut, normal):

        self.strut = strut
        self.normal = normal

        self._coordinates = None

    @property
    def coordinates(self):

        if self._coordinates is None:

            self._coordinates = (self.normal[0], self.normal[1], self.normal[2],
                                 np.dot(self.normal, self.strut))

        return self._coordinates

    def intersect(self, other):

        # compute direction of intersection line, which is the cross of
        # the normal vectors
        new_dir = np.cross(self.normal, other.normal)

        # if the new_dir is zero, then the planes are parallel
        if not np.any(new_dir):
            return None

        # compute new strut
        (a1, b1, _, d1) = self.coordinates
        (a2, b2, _, d2) = other.coordinates

        x = np.linalg.solve([[a1, b1], [a2, b2]],
                             [d1, d2])

        new_strut = (x[0], x[1], 0)

        return Line(new_strut, new_dir)


def plane_from_vec(vec1, vec2, vec3):
    """
    Computes new Plane from three vectors
    """

    # Set vec2 to be the strut, compute cross between vec2 - vec1 and vec2 - vec3 for normal
    a = np.array(vec1)
    b = np.array(vec2)
    c = np.array(vec3)

    normal = np.cross(b - a, b - c)
    return Plane(vec2, normal)
