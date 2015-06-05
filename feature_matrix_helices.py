__author__ = 'lukas'
import os
import re

from pdb.constants import  RE_HB, RE_PDB
from pdb.extract import get_secondary_structure_annotation
from learn.HydrogenBondPattern import encode_file_potential



# files with hydrogen bonds
#hb_dir = 'material/features/hydrogen_bonds'
#pdb_dir = 'material/training'

def hydrogen_bonds(hb_dir, pdb_dir):

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
    feat_pdb_paths = map(lambda x: pdb_dir + os.path.sep +x,
                         pdb_files)

    # examine all pdb files
    for (pdb, hb) in zip(feat_pdb_paths, feat_hb_paths):

        # extract helix positions
        helices = get_secondary_structure_annotation(pdb)[0]

        # create sets for alpha and 310 helices (we will ignore pi helices)
        helices_alpha = set(helices[1])
        helices_310 = set(helices[5])

        # take all entities of this pdb file into consideration
        for (w, training_point) in encode_file_potential(hb, 5):

            set_w = set(w)

            # check whether w is part of alpha helix
            if set_w <= helices_alpha:
                y = 1

            # check the same for 310 helix
            elif set_w <= helices_310:
                y = 3

            # w does not seem to belong to a helix
            else:
                 y = 0

            X.append(training_point)
            Y.append(y)

    return X, Y
