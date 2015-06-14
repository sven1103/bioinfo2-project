"""
This file is used to produce some neat evaluation results of the
ML predictor.
"""
import sys
import argparse
import os
import cPickle

from src.util import absolute_file_paths
from src.evaluation import cross_validate_split, cv_annotator
from src.conf import conf


def main(args):

    parser = argparse.ArgumentParser()
    parser.add_argument('MODE', type=str, default='CV')
    parser.add_argument('PDB_FILES', type=str)
    parser.add_argument('-f', '--fold', type=int)
    parser.add_argument('-r', type=str, default='hs')

    args = parser.parse_args(args[1:])

    conf.set_dir('.')

    pdb_files = list(absolute_file_paths(args.PDB_FILES))

    if args.MODE == 'CV':

        # split data into training and test data according to the
        # specified fold
        folds = cross_validate_split(pdb_files, args.fold)

    elif args.MODE == 'LOO':

        folds = cross_validate_split(pdb_files, len(pdb_files))

    else:
        sys.stderr.write('Invalid Validation Method' + os.linesep)
        sys.exit(3)

    accs = []

    # compue all folds
    for acc in cv_annotator(folds):
        accs.append(acc)

    with open(conf.eval_dir + os.path.sep + 'ACC', 'w') as f:
        cPickle.dump(accs, f)

if __name__ == '__main__':
    main(sys.argv)
