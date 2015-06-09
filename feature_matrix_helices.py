__author__ = 'lukas'

import os
import re

from Bio.PDB import PDBParser
from sklearn.cross_validation import  cross_val_score

from sklearn.tree import DecisionTreeClassifier

from pdb.extract import get_secondary_structure_annotation, get_id
from learn.WindowExtractor import WindowExtractor
from learn.Annotator import Annotator
from pdb.constants import RE_PDB, RE_HB
from learn.features.HydrogenBondPatternFile import get_hydrogen_bond_pattern_file
from learn.features.ChouFasmanHelix import ChouFasmanHelix




# files with hydrogen bonds
hb_dir = 'material/features/hydrogen_bonds'
pdb_dir = 'material/training'


class FeatureContext(object):

    def __init__(self, pdb_files):

        self.pdb_files = pdb_files

    @staticmethod
    def _update(pdb_path, features):
        for feature in features:
            feature.set_context(pdb_path)

    def construct_matrix(self, features, annotator, window_size):

        X = []
        Y = []

        for pdb in self.pdb_files:

            # set all features to the the current PDB contest
            FeatureContext._update(pdb, features)

            # set up Window Extractor for current PDB file
            with open(pdb, 'r') as f:

                struc = PDBParser().get_structure(get_id(pdb), f)

                helix_aa, sheet_aa = get_secondary_structure_annotation(pdb)

                we = WindowExtractor(struc, window_size, features)

                for (positions, training_point) in we.entities():

                    X.append(training_point)
                    Y.append(annotator(positions, helix_aa, sheet_aa))

        return X, Y

pdb_files = map(lambda y: pdb_dir + os.path.sep + y,
                filter(lambda x: re.match(RE_PDB, x) is not None,
                       os.listdir(pdb_dir)))

hb_files = map(lambda y: hb_dir + os.path.sep + y,
               filter(lambda x: re.match(RE_HB, x) is not None,
                      os.listdir(hb_dir)))


fc = FeatureContext(pdb_files)

hydrogen_bonds = get_hydrogen_bond_pattern_file(hb_files)()
chou_fasman = ChouFasmanHelix()

# only interested in Helices
annotator = Annotator([0, 1, 2, 3])

X, Y = fc.construct_matrix([hydrogen_bonds], annotator, 5)

print X

clf = DecisionTreeClassifier(max_depth=3)

print cross_val_score(clf, X, Y, cv=4)
