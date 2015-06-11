__author__ = 'lukas'

import os
import re

from Bio.PDB import PDBParser

from sklearn.cross_validation import  cross_val_score

from sklearn.ensemble import RandomForestClassifier

from pdb.extract import get_secondary_structure_annotation, get_id
from learn.WindowExtractor import WindowExtractor
from learn.ClassAssigner import ClassAssigner
from pdb.constants import RE_PDB, RE_HB
from learn.features.HydrogenBondPatternFile import get_hydrogen_bond_pattern_file
from learn.features.ChouFasmanHelix import ChouFasmanHelix
from learn.features.BackboneTorsionAngles import BackboneTorsionAngles






# files with hydrogen bonds
hb_dir = 'material/features/hydrogen_bonds'
pdb_dir = '/home/lukas/Dropbox/BI2_project/material/training/'


class FeatureContext(object):
    """
    Allows the construction of feature matrices by defining features
    used for encoding entities.
    """

    def __init__(self, pdb_files):

        self.pdb_files = pdb_files

    @staticmethod
    def _update(pdb_path, features):
        """
        Sets all features into the new context of the PDB file.

        :param pdb_path: Path to the PDB files that defines the context
        :param features: List of features as defined by the "features" package
        """
        for feature in features:
            feature.set_context(pdb_path)

    def construct_matrix(self, features, annotator, window_size):
        """
        Constructs feature matrix using the features and
        wanted target annotation.

        :param features: List of features as defined by the
        feature module.
        :param annotator:
        :param window_size: Window size to be used to encode the features with.
        :return: X,Y, where X encodes an entity per row with columns
        representing the features to be used.
        """

        X = []
        Y = []

        for pdb in self.pdb_files:

            # set all features to the the current PDB contest
            FeatureContext._update(pdb, features)

            # set up Window Extractor for current PDB file
            with open(pdb, 'r') as f:

                struc = PDBParser().get_structure(get_id(pdb), f)

                helix_aa, strand_aa, sheet_aa =\
                    get_secondary_structure_annotation(pdb)

                print sheet_aa

                we = WindowExtractor(struc, window_size, features)

                for (positions, training_point) in we.entities():

                    X.append(training_point)
                    Y.append(annotator(positions, helix_aa, strand_aa))

        return X, Y


def main():

    pdb_files = map(lambda y: pdb_dir + os.path.sep + y,
                    filter(lambda x: re.match(RE_PDB, x) is not None,
                           os.listdir(pdb_dir)))

    hb_files = map(lambda y: hb_dir + os.path.sep + y,
                   filter(lambda x: re.match(RE_HB, x) is not None,
                          os.listdir(hb_dir)))

    fc = FeatureContext(pdb_files)

    hydrogen_bonds = get_hydrogen_bond_pattern_file(hb_files)()
    chou_fasman = ChouFasmanHelix()
    torsion_angles = BackboneTorsionAngles()

    # only interested in Helices
    annotator = ClassAssigner([0, 1, 2, 3])

    X, Y = fc.construct_matrix([torsion_angles],
                               annotator, 5)

    clf = RandomForestClassifier(n_estimators=10)

    print cross_val_score(clf, X, Y, cv=4)


if __name__ == '__main__':
    main()
