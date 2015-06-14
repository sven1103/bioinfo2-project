from __future__ import division
from itertools import combinations
import os
import re
import cPickle

import numpy as np
from Bio.PDB import PDBParser

import src.pdb.constants as co
from src.learn.dssp import potential
from src.pdb.hydrogen import validate
from WindowFeature import WindowFeature
from SheetFeature import SheetFeature
from src.util import get_pos, get_id
from src.pdb.extract import get_amino_acids
from src.conf import conf


class HydrogenBondPattern(WindowFeature, SheetFeature):

    def __init__(self, min_seq_distance, mode='potential'):
        super(HydrogenBondPattern, self).__init__()

        # Parameters of Hydrogen Bond Patterns
        self.min_seq_distance = min_seq_distance
        self.mode = mode

        # Context of Hydrogen Bond Pattern
        self.pdb_path = None
        self.potential_map = None

    def _compute_hydrogen_bonds(self, entity):
        """
        Computes all Hydrogen Bonds that are found in the Sequence
        of Amino Acids called `entity`.

        :param entity: Iterable yielding Amino Acids
        :return: Generator of all Hydrogen Bonds found in `entity`
        """

        for (aa1, aa2) in combinations(entity, 2):

            # do not consider this pair if the number of atoms of the
            # residues is not sufficient
            if not (validate(aa1) and validate(aa2)):
                continue

            # stores both potentials between aa1 and aa2
            potentials = []

            segid1 = get_pos(aa1)
            segid2 = get_pos(aa2)

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
                yield (aa1, aa2, potentials[0], potentials[1])

    def _compute_context(self):
        """
        Computes .hb file completely from scratch if it is not
        present in the .hb directory,
        """
        # get PDBID
        pdbid = get_id(self.pdb_path)

        # initialize potential map
        self.potential_map = {}

        with open(self.pdb_path, 'r') as f:
            struc = PDBParser().get_structure(pdbid, f)

        for bond in self._compute_hydrogen_bonds(get_amino_acids(struc)):
            self.potential_map[bond[:2]] = bond[2:]

        self._compute_plain_from_map()

    def _compute_plain_from_map(self):
        """
        Writes plain hb file from dictionary
        """
        # get PDBID
        pdbid = get_id(self.pdb_path)

        # write hydrogen bonds to the plain file
        with open(conf.temp_dir + os.path.sep + pdbid + '.hb', 'w') as f:

            for pos in self.potential_map:
                a = str(pos[0])
                b = str(pos[1])
                c = str(self.potential_map[pos][0])
                d = str(self.potential_map[pos][1])
                e = [a, b, c, d]

                f.write(' '.join(e) + os.linesep)

    def _compute_map_from_plain(self, hb_path):

        # TODO Documentation
        """
        :return:
        """
        pdbid = get_id(hb_path)

        with open(hb_path, 'r') as f:

            self.potential_map = {(int(splt[0]), int(splt[1])): (float(splt[2]), float(splt[3]))
                                  for splt in (re.split(co.RE_WHITESPACE, line)
                                               for line in f)}

        with open(conf.temp_dir + os.path.sep + pdbid + '.hbc', 'w') as f:
            cPickle.dump(self.potential_map, f)


####################################################################
#####################################################################

    def tell_context(self, pdb_path):

        self.pdb_path = pdb_path

        related_files = conf.hb_dict[get_id(pdb_path)]

        hb_path = None
        hbc_path = None

        for related_file in related_files:

            if related_file.endswith('.hb'):
                hb_path = related_file

            elif related_file.endswith('.hbc'):
                hbc_path = related_file

                with open(hbc_path, 'r') as f:
                    self.potential_map = cPickle.load(f)

        if hb_path is None and hbc_path is None:
            self._compute_context()

        elif hb_path is not None:
            self._compute_map_from_plain(hb_path)

        elif hbc_path is not None:
            self._compute_plain_from_map()

    def _encode_potential(self, entity):
        """
        Extracts Hydrogen Bond Pattern from List of AA `entity`

        :param entity: List of Residue objects (which actually should be
        amino acids) that are consecutive within a protein and should
        be encoded
        :return: List of amino acid pairs that form a bond within the entity
        """
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

    def _encode_pairs(self, entity):
        """

        :param entity:
        :return:
        """
        res = []

        positions = map(lambda x: x.get_id()[1],  entity)

        # consider all possible combinations of w
        for (left, right) in combinations(positions, 2):

            # append if there is in fact an hb annotation for left, right
            if (left, right) in self.potential_map:
                res.append((left, right))
        return res

    def encode(self, entity):

        if self.mode == 'potential':
            return self._encode_potential(entity)

        elif self.mode == 'pairs':
            return self._encode_pairs(entity)
        else:
            raise ValueError("Unrecognized mode to encode Hydrogen Bonds.")

    def encode_sheet(self, strand1, strand2):

        pots = [[], [], [], []]

        for (left, right) in [(c, d)
                              for c in strand1
                              for d in strand2]:

            if (left, right) in self.potential_map:

                potentials = self.potential_map[(left, right)]
                pots[0].append(potentials[0])
                pots[1].append(potentials[1])

            if (right, left) in self.potential_map:

                potentials = self.potential_map[(right, left)]
                pots[2].append(potentials[0])
                pots[3].append(potentials[1])

        for pot in pots:
            if not pot:
                pot.append(0)

        return map(np.mean, pots)
