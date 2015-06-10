__author__ = 'sven'

import os
import re

from Bio.PDB import PDBParser
from sklearn.ensemble import RandomForestClassifier

from sklearn.cross_validation import cross_val_score

from pdb.constants import RE_PDB
from old_stuff import WindowExtractorSheets

pdb_dir_sven = "/home/sven/Git/bioinformatics2/assignment_2/pdb"
pdb_dir_lukas = "/home/lukas/Dropbox/BI2_project/material/training/"

def torsion_angles(pdb_dir):
    """
    Uses the WindowExtractorSheets class to extract torsion angle triplets.
    It will compute the features and the encoding for the classification.
    :param pdb_dir: The pdb directory
    :return: feature matrix and target vector
    """

    # Init the feature matrix and the target vector
    X = []
    y = []

    # the pdb-files sorted
    pdb_files = sorted(filter(lambda x: re.match(RE_PDB, x) is not None,
                              os.listdir(pdb_dir)))
    # make absolute paths
    pdb_paths = map(lambda x: pdb_dir + os.path.sep + x, pdb_files)

    # iterate over all pdb-files
    for pdb_file in pdb_paths:
        print "Extracting information from: ", pdb_file
        # make a PDB structure object from the file
        struc = PDBParser().get_structure("", pdb_file)
        # create a WindowExtractionSheets object.
        # It will extract the sheet relevant features from the pdb-file
        # and provide a feature matrix as well as the target vector.
        window_extractor = WindowExtractorSheets.WindowExtractorSheets(struc,
                                                                       pdb_file)
        # compute the features
        window_extractor.compute_features()
        # build the feature matrix and target vector
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

    # makes a nice feature matrix and a target vector
    for feature in feature_set:
        row = []
        for torsion_pair in feature[0]:
            row.append(torsion_pair[0])
            row.append(torsion_pair[1])
        feature_matrix.append(row)
        target_vector.append(feature[1])

    return feature_matrix, target_vector


def main():
    X, y = torsion_angles(pdb_dir_lukas)

    clf = RandomForestClassifier(n_estimators=10)

    print cross_val_score(clf, X, y, cv=4)

main()
