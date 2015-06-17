"""
This file is used to produce some neat evaluation results of the
ML predictor. Currently, it will calculate ACCs and SOV values.
"""
import sys
import argparse
import os

import numpy as np

from src.util import absolute_file_paths
from src.evaluation import cross_validate_split, cv_annotator
from src.conf import conf


def main(args):

    parser = argparse.ArgumentParser()
    parser.add_argument('MODE', type=str, default='CV')
    parser.add_argument('PDB_FILES', type=str)
    parser.add_argument('-f', '--fold', type=int)
    parser.add_argument('-m', type=str, default='Q3')

    args = parser.parse_args(args[1:])

    # set working directory of central configuration
    conf.set_dir('.')

    pdb_files = list(absolute_file_paths(args.PDB_FILES))

    # determine how many folds should be used.
    if args.MODE == 'CV':

        # split data into training and test data according to the
        # specified fold
        folds = cross_validate_split(pdb_files, args.fold)

    elif args.MODE == 'LOO':

        folds = cross_validate_split(pdb_files, len(pdb_files))

    else:
        sys.stderr.write('Invalid Validation Method' + os.linesep)
        sys.exit(3)

    # keep track of ACCs and SOVs
    accs = []
    sovs = []

    # compute all folds
    for acc, sov, in cv_annotator(folds, mode=args.m):
        accs.append(acc)
        sovs.append(sov)

    # write all results in files
    with open(conf.eval_dir + os.path.sep + 'accuracies', 'w') as f:
        for acc in accs:
            f.write(str(acc) + os.linesep)

        f.write('Mean: ' + str(np.mean(accs)) + os.linesep)

    with open(conf.eval_dir + os.path.sep + 'SOVs', 'w') as f:
        for sov in sovs:
            f.write(str(sov) + os.linesep)

        f.write('Mean: ' + str(np.mean(sovs)) + os.linesep)

if __name__ == '__main__':
    main(sys.argv)
