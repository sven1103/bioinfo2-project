from __future__ import division
from itertools import combinations

import numpy as np

from pdb.extract import get_amino_acids
import pdb.constants as co


def h_bonds(struc,
            position=None,
            min_seq_distance=2):

    for (aa1, aa2) in combinations(get_amino_acids(struc), 2):

        seqid1 = aa1.get_id()[1]
        seqid2 = aa2.get_id()[1]

        # continue if aa1 and aa2 lie out of the specified position
        if position is not None and \
                not (seqid1 in position and seqid2 in position):

            continue

        # follow the minimal distance criterion of amino acids
        if np.abs(seqid1 - seqid2) < min_seq_distance:
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

        # compute relevant distances between the atoms
        r_ON = np.linalg.norm(aa1_o_carboxyl - aa2_nitrogen)
        r_CH = np.linalg.norm(aa1_c_carboxyl - aa2_hydrogen)
        r_OH = np.linalg.norm(aa1_o_carboxyl - aa2_hydrogen)
        r_CN = np.linalg.norm(aa1_c_carboxyl - aa2_nitrogen)

        # we report a Hydrogen bond between aa1 and aa2 if the
        # potential falls under the energy cutoff
        pot = potential(r_ON, r_CH, r_OH, r_CN)

        if pot < co.HBOND_THRESHOLD:
            yield (aa1, aa2)


def potential(r_ON, r_CH, r_OH, r_CN):

    return co.Q1 * co.Q2 * ((1/r_ON) + (1/r_CH) - (1/r_OH) - (1/r_CN)) *\
        co.DIMENSIONAL_FACTOR
