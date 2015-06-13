"""
Provides utility functions that are quite general and do not fit in
any other module.
"""
from __future__ import division
from itertools import islice
from collections import defaultdict
import getopt
import re
import sys
import os

import numpy as np

from src.pdb.constants import RE_PDBID


def window(seq, n=2):
    """Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ..."""
    it = iter(seq)
    result = tuple(islice(it, n))

    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def angle(vec1, vec2, deg=True):
    """
    Computes angle between vec1 and vec2. Vectors must have the same
    dimension.

    :param vec1:  First vector.
    :param vec2:  Second vector.
    :param deg: Whether to return angle in degrees instead of radians.
    :return: Angle between both vectors
    """
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    res = np.arccos((np.dot(vec1, vec2)) / (norm1 * norm2))

    # An "expected type Union[ndarray, Iterable] might be issued here.
    # This is not a problem, as some IDEs (e.g. PyCharm) apparently do
    # not understand entire Numpy.
    return res if not deg else np.degrees(res)

def reformat_list(xs):
    """
    Reformats a list with lists of tuples into a single list of tuples
    :param xs: a list of lists with tuples
    :return: a single list of tuples
    """
    return [tuples for element in xs for tuples in element if element]


def absolute_file_paths(directory):

    for dirpath, _, filenames in os.walk(directory):

        for f in filenames:
            yield os.path.abspath(os.path.join(dirpath, f))


def pdb_map(directory):
    """
    Constructs a Dictionary from PDBID to list of
    all paths that belong to this PDBID.

    :param files: List of file paths to PDB files.
    :return: Dictionary, which maps PDBIDs to all files
    that belong to this PDBID.
    """
    res = defaultdict(list)
    for file_path in absolute_file_paths(directory):
        res[get_id(file_path)].append(file_path)

    return res

def get_pos(residue):
    """
    Returns a Position of a Residue object from the BioPython module.

    :param residue: Residue object from BioPython.
    :return: Position of Residue in AA Sequence.
    """
    return residue.get_id()[1]


def get_id(path):
    """
    Extract PDBIB from filepath.

    :param path: Path to PDB file.
    :return: PDBIDs
    """
    matches = re.findall(RE_PDBID, path)

    # if we encounter a PDBID return it, otherwise implicitely return None
    if matches:
        return matches[-1].split('.')[0]


# TODO I do not think that we actually need these functions

def printCommandLineHelp():
    print "\nBioinformatics II Project - Secondary Structure Prediction\n\
Version: 0.1\n\
Authors: Lukas Zimmermann and Sven Fillinger\n\
Required Python Version: 2.7\n\
Supported format: PDB\n\n\
Syntax: python2.7 annotate.py [options] \n\n\
Commands:\n\n\
 -h \t\t get help :)\n\
 -p [path/]  You have to provide the path to the pdb files\n" \
          " -o [path/]  output path for the calculated energies and torsions"


def check_arguments(argv):
    try:
        opts, args = getopt.getopt(argv, "hp:o:")
    except getopt.GetoptError:
        print "\n***[Error]: Wrong number of arguments.\n"
        print "***[HINT]: Did you provide the path to the pdb files?"
        printCommandLineHelp()
        sys.exit()
    if not opts:
        printCommandLineHelp()
        sys.exit()
    for opt, arg in opts:
        if opt == "-h":
            printCommandLineHelp()
            sys.exit()
        elif opt == "-p":
            path_to_pdb = arg
        elif opt == "-o":
            path_output = arg
        elif not opt:
            print "ojee"
    return path_to_pdb, path_output
