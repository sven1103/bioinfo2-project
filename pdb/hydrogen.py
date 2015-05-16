"""
This file computes explicit hydrogen atoms in PDB files by taking several rules in account, which are
described in [1].


[1] E. Baker, R. Hubbard: Hydrogen bonding in globular proteins
"""
__author__ = 'lukas'
import numpy as np

from . import extract
from algebra import plane


# TODO Currently, this funtion computs CA-N-C angles, which is not its purpose,
# This was just for testing. Should generate all explicit hydrogens of the backbone
def backbone(ident, file):

    c_carboxyl = None

    for aa in extract.get_amino_acids(ident, file):

        atoms = aa.get_unpacked_list()

        # ensures that we use the c_carboxyl of the previous iteration
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
    """
    Computes the explicit hydrogen (actually only the coordinates) of the amino acid residue `aa` and
    the c_carboxyl of the previous amino acid.
    :param aa: Coordinates of explicit hydrogen
    :param c_carboxyl:  C Carboxyl Atom of the amino acid that precedes `aa` in the sequence
    :return: Coordinates of the explicit Hydrogen of that amino acid
    """

    atoms = aa.get_unpacked_list()

    # we need Backbone Nitrogen, CA, C, 0
    nitrogen = np.array(atoms[0].get_coord())
    c_alpha = np.array(atoms[1].get_coord())
    o_carboxyl = np.array(atoms[3].get_coord())

    # setup plane through nitrogen, c_carboxyl, and o_carboxyl
    plane1 = plane.plane_from_vec(o_carboxyl,
                                  c_carboxyl,
                                  nitrogen)

    # setup bisector plane of angle C_carboxyl - N - c_alpha
    plane2 = plane.Plane(nitrogen,
                         c_carboxyl - c_alpha)

    # intersect planes to get line on which H lies
    line = plane1.intersect(plane2)

    # get candidates for positions of Hydrogen via solving an quadratic formula
    # actually, all details of ths quadratic formula have been elaborated on my flipchart. Hope there is no error here
    # TODO Check for errors, this is highly sensitive stuff here
    t = nitrogen - np.array(line.strut)

    b0 = line.direction[0]
    b1 = line.direction[1]
    b2 = line.direction[2]

    p0 = np.power(b0, 2) + np.power(b1, 2) + np.power(b2, 2)
    p1 = -1 * (2 * t[0] * b0 + 2 * t[1] * b1 + 2 * t[2] * b2)
    p2 = np.power(t[0], 2) + np.power(t[1], 2) + np.power(t[2], 2) - 0.9409

    roots = np.roots((p0, p1, p2))

    # decide for root which maximizes distance to oxygen.
    candidate1 = np.array(line.strut) + roots[0] * np.array(line.direction)
    candidate2 = np.array(line.strut) + roots[1] * np.array(line.direction)

    # compute distances of candidates to o_carboxyl
    candidate1_dist = np.linalg.norm(candidate1 - o_carboxyl)
    candidate2_dist = np.linalg.norm(candidate2 - o_carboxyl)

    return candidate2 if candidate2_dist > candidate1_dist else candidate1
