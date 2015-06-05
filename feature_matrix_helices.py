__author__ = 'lukas'
import os
import re

from Bio.PDB import PDBParser
from sklearn.tree import DecisionTreeClassifier
from sklearn.cross_validation import cross_val_score

from pdb.constants import  RE_HB, RE_PDB
from pdb.extract import get_secondary_structure_annotation
from learn.HydrogenBondPattern import encode_file_potential
from learn.WindowExtractor import  WindowExtractor
from learn.chou_fasman import ChouFasmanHelix



# files with hydrogen bonds
hb_dir = 'material/features/hydrogen_bonds'
pdb_dir = 'material/training'


def hydrogen_bonds(hb_dir, pdb_dir):
    """
    Generates Feature Matrix and class labels for hydrogen bon

    :param hb_dir:
    :param pdb_dir:
    :return:
    """

    X = []
    Y = []

    # sort files to ensure that they are equally long
    hb_files = sorted(filter(lambda x: re.match(RE_HB, x) is not None,
                             os.listdir(hb_dir)))
    pdb_files = sorted(filter(lambda x: re.match(RE_PDB, x) is not None,
                              os.listdir(pdb_dir)))

    # assemble paths
    feat_hb_paths = map(lambda x: hb_dir + os.path.sep + x,
                        hb_files)
    feat_pdb_paths = map(lambda x: pdb_dir + os.path.sep + x,
                         pdb_files)

    # examine all pdb files
    for (pdb, hb) in zip(feat_pdb_paths, feat_hb_paths):

        # extract helix positions
        helices = get_secondary_structure_annotation(pdb)[0]

        # take all entities of this pdb file into consideration
        for (w, training_point) in encode_file_potential(hb, pdb, 5):

            X.append(training_point)
            Y.append(assign_y(w, helices[1], helices[5]))

    return X, Y


def chou_fasman_parameters(pdb_dir):

        X = []
        Y = []

        pdb_files = sorted(filter(lambda x: re.match(RE_PDB, x) is not None,
                              os.listdir(pdb_dir)))

        feat_pdb_paths = map(lambda x: pdb_dir + os.path.sep + x,
                             pdb_files)

        for pdb in feat_pdb_paths:

            # extract helix positions
            helices = get_secondary_structure_annotation(pdb)[0]

            # parse structure
            with open(pdb, 'r') as f:
                struc = PDBParser().get_structure('', f)

            we = WindowExtractor(struc, 5, [ChouFasmanHelix()], None)

            # generate all entities
            for (w, training_point) in we.entities():

                X.append(training_point)
                Y.append(assign_y(w, helices[1], helices[5]))

        return X, Y


def assign_y(w, helices_alpha, helices_310):

    helices_alpha_set = set(helices_alpha)
    helices_310_set = set(helices_310)

    set_w = set(w)

    # check whether w is part of alpha helix
    if set_w <= helices_alpha_set:
        y = 1

    # check the same for 310 helix
    elif set_w <= helices_310_set:
        y = 3

    # w does not seem to belong to a helix
    else:
        y = 0

    return y


X1, Y1 = hydrogen_bonds(hb_dir, pdb_dir)
X2, Y2 = chou_fasman_parameters(pdb_dir)

X = [left + right for (left, right) in zip(X1,X2)]
Y = [left + right for (left, right) in zip(Y1,Y2)]
clf = DecisionTreeClassifier(max_depth=7)

print cross_val_score(clf, X, Y, cv=4)