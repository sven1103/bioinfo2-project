from __future__ import division
from random import shuffle
from math import floor

from sklearn.metrics import accuracy_score
import numpy as np

from learn import base
from learn.FeatureContext import FeatureContext
import configuration as conf


def train_test_split(xs, test=0.4):
    """
    Splits the list `xs` into training and test data.

    :param xs: List to be split into training and test data
    :param test: The fraction of test data to be considered
    :return: List of training and list of test data
    """
    n = len(xs)
    to = int(round(n * test))
    shuffle(xs)
    train = xs[0:to]
    test = xs[to:n]

    return train, test

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
    splits =[(xs[0:i*split_size] + xs[i*split_size + split_size:n],
              xs[i*split_size: i*split_size + split_size])
             for i in range(0, fold)]

    return splits


def evaluate_annotator(folds):
    """
    Evaluates annotation in a CV manner by training predictors on the
    test set and annotating all remaining PDB files with their
    secondary structure and evaluating on the remaining PDB files.

    :return: List of Accuracies
    """
    res = []
    # consider each fold of the CV
    for (train, test) in folds:

        # Fit predictors on training data
        fc = FeatureContext(train)
        X_train, Y_train = fc.construct_matrix(conf.helix_features,
                                              conf.helix_assigner,
                                              conf.helix_window_size)
        conf.helix_predictor.fit(X_train, Y_train)

        X_train, Y_train = fc.construct_matrix(conf.strand_features,
                                               conf.strand_assigner,
                                                conf.strand_window_size)
        conf.strand_predictor.fit(X_train, Y_train)

        # evaluate all PDB files of the test set
        test_accs = []

        for test_pdb in test:
            (map_true, map_predict) = base.annotate(test_pdb)

            # skip PDB file if we cannot annotate anything
            if map_true is None:
                continue

            test_accs.append(accuracy_score(map_true.values(),
                                            map_predict.values()))
        res.append(np.mean(test_accs))

    return res


print evaluate_annotator(cross_validate_split(conf.pdb_files, fold=4))

