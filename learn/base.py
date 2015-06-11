from __future__ import division
from collections import defaultdict

from Bio.PDB import PDBParser

import configuration.variables
from learn.features.Identity import Identity
from util import window, get_id
from learn.target_encoding import TARGET_CODES
from learn.WindowExtractor import WindowExtractor
from learn.FeatureContext import FeatureContext
import configuration as conf


def correct_helices(Y):
    """
    Tries to correct the predictions of helices, which are obviously incorrect

    :param Y: Vector of Helix types
    :return: Corrected class vector Y
    """

    for (index, (left, middle, right)) in enumerate(window(Y, 3)):

        # if the left frame is a coil and the right frame is also a coil,
        # then we we assume that the middle frame will also be a coil
        if left == TARGET_CODES['Coil'] and right == TARGET_CODES['Coil']:
            Y[index + 1] = TARGET_CODES['Coil']

        # if the left frame is a Alpha Helix and the right frame is a
        # Alpha Helix, then we assume the the middle frame also belongs
        # to a helix
        elif left == TARGET_CODES['Alpha-Helix'] \
                and right == TARGET_CODES['Alpha-Helix']:

            Y[index + 1] = TARGET_CODES['Alpha-Helix']

    return Y


def expand(struc, Y, window_size):
    """
    Annotates each Amino Acid of pdb_file with class label
    using the specified window size

    :param struc:  BioPython Structure object whose Amino Acids should be
    annotated.
    :param Y: Vector of class labels of windows
    :param window_size: Window size of windows
    :return: Class Labeling of each amino acid of pdb_file
    """

    # maps AA positions to respective class
    aa_map = defaultdict(lambda: TARGET_CODES['Coil'])

    we = WindowExtractor(struc, window_size, [Identity()])

    for (index, (positions, _)) in enumerate(we.entities()):

        # for all residues in the entity
        for pos in positions:

            # if pos is still Coil, we assume class that we see
            aa_map[pos] = Y[index]

    return aa_map


def merge_prediction(map_helix, map_strands):

    res = defaultdict(int)

    for pos in map_strands:
        res[pos] = map_strands[pos]

    for pos in map_helix:
        res[pos] = map_helix[pos]

    return res


def annotate(pdb_file):

    # get feature Matrix and true classes of PDB file
    fc = FeatureContext([pdb_file])

    # get feature Matrix for Helices of PDB file
    XHelix, YHelix = fc.construct_matrix(conf.helix_features,
                                         conf.helix_assigner,
                                         configuration.variables.helix_window_size)

    XStrand, YStrand = fc.construct_matrix(conf.strand_features,
                                           conf.strand_assigner,
                                           configuration.variables.strand_window_size)

    # if we cannot extract features from this PDB file, return None
    if not XHelix:
        return None, None

    with open(pdb_file, 'r') as f:
        struc = PDBParser().get_structure(get_id(pdb_file), f)

    # predict Helix position
    pred_helix = conf.helix_predictor.predict(XHelix)
    pred_helix = correct_helices(pred_helix)

    # predict Strand Positions
    pred_strand = conf.strand_predictor.predict(XStrand)

    # expand true and predicted helix and strand annotations
    pred_helix_expand = expand(struc, pred_helix,
                               configuration.variables.helix_window_size)
    pred_strand_expand = expand(struc, pred_strand,
                                configuration.variables.strand_window_size)

    # TODO Here we have to annotate Sheets correctly with their orientation
    # first, get all encompassed Hydrogen bonds
    #hydrogen_bonds = StrucExtractor(struc).get_hydrogen_bonds()

    true_helix_expand = expand(struc, YHelix,
                               configuration.variables.helix_window_size)
    true_strand_expand = expand(struc, YStrand,
                                configuration.variables.strand_window_size)

    # merge separate predictions of helices and sheets
    pred_expand = merge_prediction(pred_helix_expand, pred_strand_expand)
    true_expand = merge_prediction(true_helix_expand, true_strand_expand)

    # return helix annotation, together with its accuracy
    return (true_expand, pred_expand)
