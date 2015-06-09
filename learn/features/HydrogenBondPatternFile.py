__author__ = 'lukas'
import re
from itertools import combinations

import pdb.constants as co
from pdb.extract import get_id
from Feature import Feature


def get_hydrogen_bond_pattern_file(hb_files):

    def _construct_potential_map(hb_path):

        with open(hb_path, 'r') as f:

            # incomprehensible Code
            hydrogen_bonds = ((int(splt[0]), int(splt[1]),
                               float(splt[2]), float(splt[3]))
                              for splt in [re.split(co.RE_WHITESPACE, line)
                                           for line in f])

        # assemble potential map... somehow
        potential_map = {(aa1_pos, aa2_pos): (pot1, pot2)
                         for (aa1_pos, aa2_pos, pot1, pot2)
                         in hydrogen_bonds}

        return potential_map

    # construct map of all potential maps
    pdb_potentials = {get_id(hb_path): _construct_potential_map(hb_path)
                      for hb_path in hb_files}

    class HydrogenBondPatternFile(Feature):

        def __init__(self):
            super(HydrogenBondPatternFile, self).__init__()
            self.potential_map = None

        def set_context(self, pdb_path):
            self.potential_map = pdb_potentials[get_id(pdb_path)]

        def encode(self, entity):

            res = []

            positions = map(lambda x: x.get_id()[1],  entity)

            # consider all possible combinations of w
            for (left, right) in combinations(positions, 2):

                # append if there is in fact an hb annotation for left, right
                if (left, right) in self.potential_map:

                    # fetch potentials
                    potentials = self.potential_map[(left, right)]
                    res.extend([potentials[0], potentials[1]])
                else:
                    res.extend([0, 0])

            return res

    return HydrogenBondPatternFile
