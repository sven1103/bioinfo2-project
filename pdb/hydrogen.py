"""
This file computes explicit hydrogen atoms in PDB files by taking several rules in account, which are
described in [1].


[1] E. Baker, R. Hubbard: Hydrogen bonding in globular proteins
"""
__author__ = 'lukas'
import numpy as np

from . import extract
from algebra import plane


def backbone(ident, file):

    c_carboxyl = None

    for aa in extract.get_amino_acids(ident, file):

        atoms = aa.get_unpacked_list()

        if c_carboxyl is None:
            c_carboxyl = np.array(atoms[2].get_coord())
            continue

        # compute angle
        n = np.array(atoms[0].get_coord())
        c_a = np.array(atoms[1].get_coord())
        h = np.array(_compute_hydrogen(aa, c_carboxyl))

        v1 = c_a - n
        v2 = h - n

        angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

        yield (angle / (2 * np.pi)) * 360

        # set new c_carboxyl
        c_carboxyl = np.array(atoms[2].get_coord())


def _compute_hydrogen(aa, c_carboxyl):

    atoms = aa.get_unpacked_list()

    # we need Backbone Nitrogen, CA, C, 0
    nitrogen = np.array(atoms[0].get_coord())
    c_alpha = np.array(atoms[1].get_coord())
    o_carboxyl = np.array(atoms[3].get_coord())

    # setup plane through nitrogen, c_carboxyl, and o_carboxyl
    plane1 = plane.plane_from_vec(o_carboxyl,
                                  c_carboxyl,
                                  nitrogen)

    #TODO This approach does not work yet. Should use the angle of 119 deg in CA-N-H
    # System instead of bisector.

    # setup bisector plane of angle C_carboxyl - N - c_alpha
    plane2 = plane.Plane(nitrogen,
                         c_carboxyl - c_alpha)

    # intersect planes to get line on which H lies
    line = plane1.intersect(plane2)

    # get candidates for positions of Hydrogen via solving an quadratic formula

    t = nitrogen - np.array(line.strut)

    b0 = line.direction[0]
    b1 = line.direction[1]
    b2 = line.direction[2]

    p0 = np.power(b0, 2) + np.power(b1, 2) + np.power(b2, 2)
    p1 = -1 * (2 * t[0] * b0 + 2 * t[1] * b1 + 2 * t[2] * b2)
    p2 = np.power(t[0], 2) + np.power(t[1], 2) + np.power(t[2], 2) - 0.81

    roots = np.roots((p0, p1, p2))

    # decide for root which maximizes distance to oxygen.
    candidate1 = np.array(line.strut) + roots[0] * np.array(line.direction)
    candidate2 = np.array(line.strut) + roots[1] * np.array(line.direction)

    # compute
    candidate1_dist = np.linalg.norm(candidate1 - o_carboxyl)
    candidate2_dist = np.linalg.norm(candidate2 - o_carboxyl)

    return candidate1 if candidate1_dist < candidate2_dist else candidate2
