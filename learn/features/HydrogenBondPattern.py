from __future__ import division
from itertools import combinations

import numpy as np

import pdb.constants as co
from learn.dssp import potential
from pdb.hydrogen import validate
from Feature import Feature

from Bio import PDB
from pdb.extract import get_amino_acids
from pdb.hydrogen import StructureHydrogenAdder

class HydrogenBondPattern(Feature):

    def __init__(self, min_seq_distance):
        super(HydrogenBondPattern, self).__init__()
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
            if not (validate(aa1) and validate(aa2)):
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


if __name__ == "__main__":
    #pdb_file = "/home/sven/Git/bioinformatics2/assignment_2/pdb/1SMC.pdb"
    pdb_file = "/home/fillinger/git/bioinformatics2/assignment_2/pdb/1q4k.pdb"
    struct = PDB.PDBParser().get_structure("test", pdb_file)
    residues = get_amino_acids(struct)
    print "uhu"
    sha = StructureHydrogenAdder(struct)
    sha.supplement()
    hbonds = HydrogenBondPattern(2)
    for (aa1, aa2, pot1, pot2) in hbonds.encode(get_amino_acids(sha.struc)):

        line = str(aa1.get_id()[1]) + ' ' + str(aa2.get_id()[1]) + ' ' + str(pot1) + ' ' + str(pot2)
        print line