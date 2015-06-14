from __future__ import division
from math import floor

from sklearn.metrics import accuracy_score
import numpy as np

from src.learn import base
from learn.FeatureContext import FeatureContext
from src.conf import conf
from feature_list import helix_features, strand_features
from src.learn.target_encoding import Q3_MAPPING, HELIX_MAPPING, STRAND_MAPPING


def train_test_split(xs, test=0.4):
    """
    Splits the list `xs` into training and test data.

    :param xs: List to be split into training and test data
    :param test: The fraction of test data to be considered
    :return: List of training and list of test data
    """
    n = len(xs)
    to = int(round(n * test))
    return xs[0:to], xs[to:n]

def cross_validate_split(xs, fold=4):
    """
    Generates CV splits of the list `xs`.

    :param xs: The list to be iteratively split intro training and test data
    :param fold: Fold of Cross-Validation
    :return: List of tuples, where the ith element encompasses the training
    set in the first and the test set in the second component for split i.
    """
    n = len(xs)
    split_size = int(floor(n / fold))

    # generate all CV splits of the list xs
    splits = [(xs[0:i*split_size] + xs[i*split_size + split_size:n],
               xs[i*split_size: i*split_size + split_size])
              for i in range(0, fold)]

    return splits


def eval_pdb(test_pdb):
    """
    Evaluates the PDB file given by its path when predictor is trained
    on the training data

    :param test_pdb:
    :return:
    """
    map_true, map_predict, true_sheets, predicted_sheets\
        = base.annotate(test_pdb)

    acc = accuracy_score(map_true.values(), map_predict.values())

    return acc

def cv_annotator(folds, mode='Q3'):
    """
    Evaluates annotation in a CV manner by training predictors on the
    test set and annotating all remaining PDB files with their
    secondary structure and evaluating on the remaining PDB files.

    :return: List of Accuracies
    """
    # consider each fold of the CV
    for (train, test) in folds:

        conf.reset_predictors()

        # Fit predictors on training data
        fc = FeatureContext(train)

        if mode == 'Q3':
            mapping = Q3_MAPPING
        elif mode == 'H':
            mapping = HELIX_MAPPING
        elif mode == 'E':
            mapping = STRAND_MAPPING

        print "Start training helices"
        x_train, y_train = fc.construct_window_matrix(helix_features,
                                                      conf.helix_assigner,
                                                      conf.helix_window_size)

        conf.helix_predictor.fit(x_train, y_train)

        print "Start training Strands"
        x_train, y_train = fc.construct_window_matrix(strand_features,
                                                      conf.strand_assigner,
                                                      conf.strand_window_size)

        conf.strand_predictor.fit(x_train, y_train)
        print "Finished Training"

        acc = []
        sov = []
        # take average score of test set
        for test_pdb in test:

            print "Test", test_pdb

            map_true, map_predict = base.annotate_no_sheet(test_pdb)

            # skip PDB file if we cannot annotate anything
            if map_true is None:
                continue

            # map down to the Q3 score
            map_true = map(lambda x: mapping[x], map_true.values())
            map_predict = map(lambda x: mapping[x], map_predict.values())

            acc.append(accuracy_score(map_true, map_predict))

        yield np.mean(acc)
