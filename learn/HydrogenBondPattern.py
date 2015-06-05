from __future__ import division
from itertools import combinations
import re
from collections import defaultdict

import numpy as np

import pdb.constants as co
from learn.dssp import potential
from pdb.hydrogen import _validate
from util import window

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

            # stores both potentials between aa1 and aa2
            potentials = []

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

            # extract atoms from both amino acids
            atoms = [aa1.get_unpacked_list(),
                     aa2.get_unpacked_list()]

            for i in range(0, len(atoms)):
                c_carboxyl = np.array(atoms[i][2].get_coord())
                o_carboxyl = np.array(atoms[i][3].get_coord())

                nitrogen = np.array(atoms[1-i][0].get_coord())
                hydrogen = None
                for atom in atoms[1-i]:
                    if atom.get_name().strip() == 'H':
                        hydrogen = np.array(atom.get_coord())

                if hydrogen is None:
                    potentials.append(0)
                    continue

                # compute relevant distances
                r_ON = np.linalg.norm(o_carboxyl - nitrogen)
                r_CH = np.linalg.norm(c_carboxyl - hydrogen)
                r_OH = np.linalg.norm(o_carboxyl - hydrogen)
                r_CN = np.linalg.norm(c_carboxyl - nitrogen)

                # compute potential
                pot = potential(r_ON, r_CH, r_OH, r_CN)

                potentials.append(pot if pot < co.HBOND_THRESHOLD else 0)

            # we return this as an result if at least one potential
            # is below the threshold , so they are not both 0

            if sum(potentials) != 0:
                res.append((aa1, aa2, potentials[0], potentials[1]))

        return res

def read_window(file_path, window_size, positions=None):
    """
    Yields all entities of length `window_size` for the specified
    file in `file_path` only for the Positions `positions` we are interested
    in.

    :param file_path:
    :param window_size:
    :param positions:
    :return:
    """
    # fetch all hydrogen bonds in the specified file
    with open(file_path, 'r') as f:

        hydrogen_bonds = []

        for line in f:

            line = line.strip()

            # split line for whitespace
            splt = re.split(co.RE_WHITESPACE, line)
            hydrogen_bonds.append((int(splt[0]), int(splt[1]),
                                   float(splt[2]), float(splt[3])))

        # extract windows from the hydrogen bonds
        # determine max position of hydrogen bonds
        max_position = max(map(lambda x: max(x[0], x[1]),
                               hydrogen_bonds))
        print max_position

        # map all position tuples to their potentials
        potential_map = defaultdict(list)

        for (aa1_pos, aa2_pos, pot1, pot2) in hydrogen_bonds:
            potential_map[(aa1_pos, aa2_pos)] = (pot1, pot2)

        # use window of size window_size
        for w in window(range(1, max_position+1), window_size):

            entity = []

            # consider all possible combinations of w
            for (left, right) in combinations(w, 2):

                # append if there is in fact an hb annotation for left, right
                if (left, right) in potential_map:

                    # fetch potentials
                    potentials = potential_map[(left, right)]
                    entity.append((left, right, potentials[0], potentials[1]))

            yield entity