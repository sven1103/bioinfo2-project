__author__ = 'sven'

import os
import re

from Bio.PDB import PDBParser

from pdb.constants import RE_PDB
from learn import WindowExtractorSheets

from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import cross_val_score

pdb_dir = "/home/sven/Git/bioinformatics2/assignment_2/pdb"


def torsion_angles(pdb_dir):

    # Init the feature matrix and the target vector
    X = []
    y = []

    pdb_files = sorted(filter(lambda x: re.match(RE_PDB, x) is not None,
                              os.listdir(pdb_dir)))

    pdb_paths = map(lambda x: pdb_dir + os.path.sep + x, pdb_files)
    print pdb_paths

    for pdb_file in pdb_paths[1:]:
        print "Extracting information from: ", pdb_file
        struc = PDBParser().get_structure("", pdb_file)
        window_extractor = WindowExtractorSheets.WindowExtractorSheets(struc, pdb_file)
        window_extractor.compute_features()
        X_temp, y_temp = feature_converter(window_extractor)
        X.extend(X_temp)
        y.extend(y_temp)

    return X, y


def feature_converter(window_extractor_obj):
    """
    Makes a feature matrix and target vector out of the feature list
    from an WindowExtractorSheet object
    :param window_extractor_obj: a WindowExtractorSheets-object
    :return: the feature matrix and the target vector
    """
    feature_set = window_extractor_obj.get_features()

    feature_matrix = []
    target_vector = []

    for feature in feature_set:
        row = []
        for torsion_pair in feature[0]:
            row.append(torsion_pair[0])
            row.append(torsion_pair[1])
        feature_matrix.append(row)
        target_vector.append(feature[1])

    return feature_matrix, target_vector


def main():
    X, y = torsion_angles(pdb_dir)

    clf = RandomForestClassifier(n_estimators=10)

    print cross_val_score(clf, X, y, cv=4)

main()
