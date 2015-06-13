"""
This file is used to create new Prediction Models that will then be loaded
into the annotate.py Script.
"""
import sys
import argparse
import os
import cPickle

from src.util import absolute_file_paths
from src.conf import conf
from feature_list import helix_features, strand_features, sheet_features
from src.learn.FeatureContext import FeatureContext

PARSER_DESC = ("Accepts a PDB file and creates Predictors based on "
               "the parameters set in configuration.py and feature_list.py")

def main(argv):

    parser = argparse.ArgumentParser(description=PARSER_DESC)
    parser.add_argument('PDB_FILES', type=str)

    args = parser.parse_args(argv[1:])
    conf.reset_predictors()
    conf.set_dir('.')

    pdb_files = absolute_file_paths(args.PDB_FILES)

    # Get Feature table of Training data
    fc = FeatureContext(pdb_files)

    # train Helix Predictor
    x_train, y_train = fc.construct_window_matrix(helix_features,
                                                 conf.helix_assigner,
                                                 conf.helix_window_size)
    conf.helix_predictor.fit(x_train, y_train)

    # train Strand Predictor
    y_train, y_train = fc.construct_window_matrix(strand_features,
                                                 conf.strand_assigner,
                                                 conf.strand_window_size)
    conf.strand_predictor.fit(x_train, y_train)

    # train Sheet Predictor
    x_train, y_train = fc.construct_sheet_matrix(sheet_features)

    conf.sheet_predictor.fit(x_train, y_train)

    with open(conf.pred_dir + os.path.sep + 'HELIX', 'w') as f:
        cPickle.dump(conf.helix_predictor, f)
    with open(conf.pred_dir + os.path.sep + 'HELIX', 'w') as f:
        cPickle.dump(conf.helix_predictor, f)
    with open(conf.pred_dir + os.path.sep + 'HELIX', 'w') as f:
        cPickle.dump(conf.helix_predictor, f)



if __name__ == '__main__':
    main(sys.argv)