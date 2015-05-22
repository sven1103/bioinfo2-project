from __future__ import division

import numpy as np

from algebra.Line import Line


class Plane(object):
    """
    Represents plane in 3D space.
    """

    def __init__(self, strut, normal):

        self.strut = np.array(strut)
        self.normal = np.array(normal)

        self._coordinates = None

    @property
    def coordinates(self):

        if self._coordinates is None:

            self._coordinates = (self.normal[0], self.normal[1],
                                 self.normal[2],
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


def from_vec(vec1, vec2, vec3):
    """
    Computes new Plane that is defined by 3 points in 3D space.

    :param vec1: First point on new Plane.
    :param vec2: Second point on new Plane.
    :param vec3: Third Point on new Plane.
    :return: Plane object containing all three argument points.
    """

    # Set vec2 to be the strut, compute cross between vec2 - vec1
    # and vec2 - vec3 for normal
    a = np.array(vec1)
    b = np.array(vec2)
    c = np.array(vec3)

    normal = np.cross(b - a, b - c)
    return Plane(vec2, normal)
