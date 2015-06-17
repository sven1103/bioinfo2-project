from __future__ import division
from math import floor

from sklearn.metrics import accuracy_score
import numpy as np

from src.learn import base
from learn.FeatureContext import FeatureContext
from src.conf import conf
from feature_list import helix_features, strand_features
from src.learn.target_encoding import Q3_MAPPING, HELIX_MAPPING, STRAND_MAPPING


# extract segments:
# extract segments:
def get_segments(positions, segment_symbol):

    segments = []

    in_segment = False

    for (index, symbol) in enumerate(positions):

        if symbol == segment_symbol and not in_segment:
            in_segment = True
            left = index
        elif symbol != segment_symbol and in_segment:
            in_segment = False
            right = index - 1
            segments.append((left, right))

    if in_segment:
        right = len(positions) - 1
        segments.append((left, right))

    return segments


def minov(seg1, seg2):

    return len(set(range(seg1[0], seg1[1] + 1)) &
               set(range(seg2[0], seg2[1] + 1)))


def maxov(seg1, seg2):

    return len(set(range(seg1[0], seg1[1] + 1)) |
               set(range(seg2[0], seg2[1] + 1)))

def seglen(seg):

    return seg[1] - seg[0] + 1


def delta(seg1, seg2):

    return min([maxov(seg1, seg2) - minov(seg1, seg2),
               minov(seg1,seg2),
               seglen(seg1) / 2,
               seglen(seg2) / 2])


def overlap(seg1, seg2):

    return minov(seg1, seg2) != 0


def generate_overlapping(segments1, segments2):

    for (seg1, seg2) in ((c, d)
                         for c in segments1
                         for d in segments2):

        if overlap(seg1, seg2):
            yield seg1, seg2




def segment_overlap(positions1, positions2, types):

    segments1_list = [get_segments(positions1, segment_symbol=t)
                      for t in types]
    segments2_list = [get_segments(positions2, segment_symbol=t)
                      for t in types]

    n = sum((sum(map(seglen, segments)) for segments in segments1_list))

    return (1/n) * sum((sum((((minov(seg1, seg2) + delta(seg1, seg2)) /
                              (maxov(seg1, seg2))) * seglen(seg1)
                             for (seg1, seg2) in generate_overlapping(segments1,
                                                                      segments2)))
                        for (segments1, segments2) in zip(segments1_list,
                                                          segments2_list)))


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
            symbols = ['H', 'E', '-']

        elif mode == 'H':
            mapping = HELIX_MAPPING
            symbols = ['H', '-']

        elif mode == 'E':
            mapping = STRAND_MAPPING
            symbols = ['E', '-']

        x_train, y_train = fc.construct_window_matrix(helix_features,
                                                      conf.helix_assigner,
                                                      conf.helix_window_size)

        conf.helix_predictor.fit(x_train, y_train)

        x_train, y_train = fc.construct_window_matrix(strand_features,
                                                      conf.strand_assigner,
                                                      conf.strand_window_size)

        conf.strand_predictor.fit(x_train, y_train)

        acc = []
        sov = []

        # take average score of test set
        for test_pdb in test:

            map_true, map_predict = base.annotate_no_sheet(test_pdb)

            # skip PDB file if we cannot annotate anything
            if map_true is None:
                continue

            # map down to the Q3 score
            map_true = map(lambda x: mapping[x], map_true.values())
            map_predict = map(lambda x: mapping[x], map_predict.values())

            acc.append(accuracy_score(map_true, map_predict))
            sov.append(segment_overlap(map_true, map_predict, symbols))

        yield np.mean(acc), np.mean(sov),
