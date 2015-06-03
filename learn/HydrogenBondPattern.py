from __future__ import division
from itertools import combinations

import numpy as np

import pdb.constants as co
from learn.dssp import potential
from pdb.hydrogen import _validate


class HydrogenBondPattern(object):

    def __init__(self, min_seq_distance):
        self.min_seq_distance = min_seq_distance

    def encode(self, entity):
        """
        Extracts Hydrogen Bond Pattern from List of AA `entity`

        :param entity: List of Residue objects (which actually should be
        amino acids) that are consecutive within a protein and should
        be encoded
        :return: List of amino acid pairs that form a bond within the entity
        """
        # generate all unique pairs of amino acids within this entity
        res = []

        for (aa1, aa2) in combinations(entity, 2):

            # do not consider this pair if the number of atoms of the
            # residues is not sufficient
            if not (_validate(aa1) and _validate(aa2)):
                continue

            segid1 = aa1.get_id()[1]
            segid2 = aa2.get_id()[1]

            # distance
            dist = np.abs(segid1 - segid2)

            # take care of the minimal sequence distance criterion
            # between aa1 and aa2
            if dist < self.min_seq_distance:
                continue

            # get relevant atom positions
            atoms1 = aa1.get_unpacked_list()
            atoms2 = aa2.get_unpacked_list()

            aa1_c_carboxyl = np.array(atoms1[2].get_coord())
            aa1_o_carboxyl = np.array(atoms1[3].get_coord())
            aa2_nitrogen = np.array(atoms2[0].get_coord())

            # find hydrogen in atoms2
            for atom in atoms2:
                if atom.get_name().strip() == 'H':
                    aa2_hydrogen = np.array(atom.get_coord())

            r_ON = np.linalg.norm(aa1_o_carboxyl - aa2_nitrogen)
            r_CH = np.linalg.norm(aa1_c_carboxyl - aa2_hydrogen)
            r_OH = np.linalg.norm(aa1_o_carboxyl - aa2_hydrogen)
            r_CN = np.linalg.norm(aa1_c_carboxyl - aa2_nitrogen)
            pot = potential(r_ON, r_CH, r_OH, r_CN)

            if pot < co.HBOND_THRESHOLD:
                res.append((aa1, aa2, pot))

        return res

def read_window(file_path, window_size, positions=None):
    """
    Returns all entities of length `window_size` for the specified
    file in `file_path` only for the Positions `positions` we are interested
    in.

    :param file_path:
    :param window_size:
    :param positions:
    :return:
    """
    pass
