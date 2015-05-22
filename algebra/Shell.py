"""
This file provides a class for a sphere in 3D, but only
representing the outer shell

"""
import numpy as np

import Plane


class Shell(object):

    def __init__(self, center, radius):

        self.center = np.array(center)
        self.radius = radius

    def element(self, theta, phi):
        """
        Returns vector that is element of this shell by specifying
        angles theta and phi

        :param theta: First angle of vector
        :param phi: Second angle of vector
        :return: Vector which is element of this shell.
        """

        a = np.sin(theta)
        b = np.cos(theta)
        c = np.sin(phi)
        d = np.cos(phi)

        angles = np.array((b * c, a * c, d))

        return self.center + (self.radius * angles)

    def intersect(self, other):
        """
        Intersect this shell with a plane, yields function
        which return 0 iff vector specified by theta and
        phi lies on the plane

        :param other:
        :return:
        """

        if not isinstance(other, Plane.Plane):
            raise AttributeError(("You cannot intersect a Shell with this"
                                 "object."))

        def equation(theta, phi):

            # define all components of the equation
            n1 = other.normal[0]
            n2 = other.normal[1]
            n3 = other.normal[2]
            x0 = self.center[0]
            y0 = self.center[1]
            z0 = self.center[2]
            p1 = other.strut[0]
            p2 = other.strut[1]
            p3 = other.strut[2]
            r = self.radius

            # setup parts of the equation
            a = n1 * x0 - p1 * n1 + n2 * y0 - p2 * n2 + n3 * z0 - p3 * n3
            b = n1 * r * np.cos(theta) * np.sin(phi)
            c = n2 * r * np.sin(theta) * np.sin(phi)
            d = n3 * r * np.cos(phi)

            return a + b + c + d

        return equation
