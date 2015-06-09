"""
This file computes explicit hydrogen atoms in PDB files by taking several rules
in account, which are described in [1].

[1] E. Baker, R. Hubbard: Hydrogen bonding in globular proteins
"""
from __future__ import division

from Bio.PDB import Atom
import numpy as np
from scipy.optimize import fmin_slsqp

import pdb.constants
from util import window, angle
import extract
from algebra import Plane


class StructureHydrogenAdder(object):

    def __init__(self, struc):

        self.struc = struc

        self._atom_max_serial_number = extract.get_atom_max_serial_number(struc)
        self._current_atom_serial_number = self._atom_max_serial_number

    def _next_atom_serial(self):
        """
        Returns next atom serial that can be used as serial number
        for a newly generated atom.

        :return: Next free atom serial number.
        """

        self._current_atom_serial_number += 1
        return self._current_atom_serial_number

    def _aa_window(self, n):
        """
        Generates all consecutive pairs of amino acids, provided
        that the aa are valid (which means that they consist of at the four
        atoms.
        """
        for (aa1, aa2) in window(extract.get_amino_acids(self.struc), n):

            if validate(aa1) and validate(aa2):
                yield (aa1, aa2)

    def supplement(self):
        """
        Adds all missing hydrogen to this structure.
        """

        for (aa1, aa2) in self._aa_window(2):

            # add hydrogen to residue, i.e. aa2
            # only compute hydrogen if there is not already a hydrogen defined
            if not _contains_hydrogen(aa2):
                aa2.add(self._compute_hydrogen(aa1, aa2))

    def backbone(self):
        """
        Generates all explicit Hydrogens of Backbone of structure
        object. Does not check whether these are already included in
        the residues.

        :return: 2-Tuple (Hydrogen, resid), where Hydrogen is an object of class
        Atom that represents one computed hydrogen. resid is the ID of
        the Entity object that contains the hydrogen.
        """

        # uses a sliding window approach to generate all consecutive
        # tuples of amino acids
        for (aa1, aa2) in self._aa_window(2):

            # compute explicit hydrogen of aa2
            aa2_hydrogen = self._compute_hydrogen(aa1, aa2)

            yield (aa2_hydrogen, aa2.get_full_id())

    def _compute_hydrogen(self, aa1, aa2):
        """
        Computes the explicit hydrogen (actually only the coordinates) of
        the amino acid residue `aa` and the c_carboxyl of the previous amino acid.
        :param aa1: Amino acid one
        :param aa2: Amino acid two
        :return: Coordinates of the explicit Hydrogen of amino acid aa2
        """
        # retrieve all atoms of considered amino acids
        atoms1 = aa1.get_unpacked_list()
        atoms2 = aa2.get_unpacked_list()

        # get coordinates of all relevant atoms
        aa1_oxygen = np.array(atoms1[3].get_coord())
        aa1_carbon = np.array(atoms1[2].get_coord())
        aa2_nitrogen = np.array(atoms2[0].get_coord())
        aa2_calpha = np.array(atoms2[1].get_coord())

        # setup plane through aa1_oxygen, aa1_carbon, and aa2_nitrogen
        plane1 = Plane.from_vec(aa1_oxygen, aa1_carbon, aa2_nitrogen)

        # setup bisector plane
        plane2 = Plane.Plane(aa2_nitrogen,
                             aa2_calpha - aa1_carbon)

        # intersect planes to get line on which H lies
        line = plane1.intersect(plane2)

        # if the planes do not intersect, we we set a simple starting
        # position for optimization
        if line is None:
            hydrogen = np.array((0, 0, 0))

        else:
            # get candidates for positions of Hydrogen
            # via solving an quadratic formula
            t = aa2_nitrogen - line.strut

            b0 = line.direction[0]
            b1 = line.direction[1]
            b2 = line.direction[2]

            p0 = np.power(b0, 2) + np.power(b1, 2) + np.power(b2, 2)
            p1 = (-1) * (2 * t[0] * b0 + 2 * t[1] * b1 + 2 * t[2] * b2)
            p2 = np.power(t[0], 2) + np.power(t[1], 2) + np.power(t[2], 2)\
                - pdb.constants.NH_DISTANCE2

            roots = np.roots((p0, p1, p2))

            # decide for root which maximizes distance to oxygen.
            candidate1 = line.strut + roots[0] * line.direction
            candidate2 = line.strut + roots[1] * line.direction

            # compute distances of candidates to oxygen of aa1
            candidate1_dist = np.linalg.norm(candidate1 - aa1_oxygen)
            candidate2_dist = np.linalg.norm(candidate2 - aa1_oxygen)

            hydrogen = candidate2 if candidate2_dist > candidate1_dist \
                else candidate1

        #
        # optimize hydrogen angle
        #
        # get coordinates of first plane
        p1_coord = plane1.coordinates

        # ensure that the hydrogen lies on the plane
        def on_plane(x, y, z):
            return p1_coord[0] * x + p1_coord[1] * y + p1_coord[2] * z\
                - p1_coord[3]

        # ensures the correct distance between nitrogen and hydrogen
        def has_distance(x, y, z):
            return np.linalg.norm(vec1(x, y, z)) - pdb.constants.NH_DISTANCE

        # consider vectors from nitrogen to hydrogen and vector
        # from calpha to nitrogen
        vec1 = lambda x, y, z: np.array((x, y, z)) - aa2_nitrogen
        vec2 = aa2_calpha - aa2_nitrogen

        # specify function that is to be minimized
        def target(x, y, z):

            return np.power(angle(vec1(x, y, z), vec2, deg=True)
                            - pdb.constants.HNCA_ANGLE_DEG, 2)

        res = fmin_slsqp(lambda q: target(q[0], q[1], q[2]),
                         x0=hydrogen,
                         eqcons=[lambda q: on_plane(q[0], q[1], q[2]),
                                 lambda q: has_distance(q[0], q[1], q[2])],
                         acc=1e-12,
                         epsilon=1e-10,
                         iprint=0,
                         iter=1000)

        # create Atom instance of computed Hydrogen and return
        return Atom.Atom('H', res, 0, 1, ' ', ' H  ',
                         self._next_atom_serial(), 'H')

def validate(aa):
    """
    Check whether amino acid aa is valid (it must contain at least
    four atoms. Otherwise, something is wrong.

    :param aa: Amino acid (residue in BioPython) to be validated.
    :return: True iff aa is valid
    """

    return len(aa.get_unpacked_list()) > 3


def _contains_hydrogen(aa):
    """
    Checks whether aa contains backbone hydrogen
    """
    for atom in aa.get_unpacked_list():

        if atom.get_name().strip() == 'H':
            return True
    return False
