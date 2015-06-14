from __future__ import division
from collections import defaultdict

from Bio.PDB import PDBParser

from src.learn.features.WindowIdentity import WindowIdentity
from src.learn.sheets import predict_sheets
from src.util import window, get_id
from src.learn.target_encoding import TARGET_CODES
from src.learn.WindowExtractor import WindowExtractor
from src.learn.FeatureContext import FeatureContext
from src.learn.features.HydrogenBondPattern import HydrogenBondPattern
from src.pdb.extract import get_amino_acids
from src.conf import conf
from feature_list import helix_features, strand_features

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

def correct_strand(Y):
    """
    Tries to correct the predictions of strands, which are obviously incorrect

    :param Y: Vector of Strand types
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
        elif left == TARGET_CODES['sss'] \
                and right == TARGET_CODES['sss']:

            Y[index + 1] = TARGET_CODES['sss']

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

    we = WindowExtractor(struc, window_size, [WindowIdentity()])

    for (index, (positions, _)) in enumerate(we.entities()):

        # for all residues in the entity
        for pos in positions:

            # if pos is still Coil, we assume class that we see
            if aa_map[pos] == TARGET_CODES['Coil']:
                aa_map[pos] = Y[index]

    return aa_map


def merge_prediction(map_helix, map_strands):

    res = defaultdict(int)

    for pos in map_strands:
        res[pos] = map_strands[pos]

    for pos in map_helix:

        if map_helix[pos] != TARGET_CODES['Coil']:

            res[pos] = map_helix[pos]

    return res


def annotate(pdb_file):

    # get feature Matrix and true classes of PDB file
    fc = FeatureContext([pdb_file])

    # get feature Matrix for Helices of PDB file
    XHelix, YHelix = fc.construct_window_matrix(helix_features,
                                                conf.helix_assigner,
                                                conf.helix_window_size)

    XStrand, YStrand = fc.construct_window_matrix(strand_features,
                                                  conf.strand_assigner,
                                                  conf.strand_window_size)

    # if we cannot extract features from this PDB file, return None
    if not XHelix:
        return None, None

    with open(pdb_file, 'r') as f:
        struc = PDBParser().get_structure(get_id(pdb_file), f)

    # predict Helix position
    helix_predictor = conf.helix_predictor
    pred_helix = helix_predictor.predict(XHelix)
    pred_helix = correct_helices(pred_helix)

    # predict Strand Positions
    strand_predictor = conf.strand_predictor
    pred_strand = strand_predictor.predict(XStrand)
    pred_strand = correct_strand(pred_strand)

    # expand predicted helix and strand annotations
    pred_helix_expand = expand(struc, pred_helix,
                               conf.helix_window_size)
    pred_strand_expand = expand(struc, pred_strand,
                                conf.strand_window_size)

    # first, get all encompassed Hydrogen bonds

    sheet_bonds = HydrogenBondPattern(2, mode='pairs')
    sheet_bonds.tell_context(pdb_file)
    hydrogen_bonds = sheet_bonds.encode(get_amino_acids(struc))

    # get true annotations of sheets
    true_sheets = fc.get_sheets()[pdb_file]

    predicted_strand_positions = [k for (k, v) in pred_strand_expand.iteritems()
                        if v != TARGET_CODES['Coil']]

    pred_sheets = predict_sheets(hydrogen_bonds, predicted_strand_positions,
                                 pdb_file)

    true_helix_expand = expand(struc, YHelix,
                               conf.helix_window_size)
    true_strand_expand = expand(struc, YStrand,
                                conf.strand_window_size)

    true_strand_positions = [k for (k, v) in true_strand_expand.iteritems()
                        if v != TARGET_CODES['Coil']]

    print true_strand_positions
    print predicted_strand_positions

    # merge separate predictions of helices and sheets
    pred_expand = merge_prediction(pred_helix_expand, pred_strand_expand)
    true_expand = merge_prediction(true_helix_expand, true_strand_expand)

    return true_expand, pred_expand, true_sheets, pred_sheets


def annotate_no_sheet(pdb_file):

    # get feature Matrix and true classes of PDB file
    fc = FeatureContext([pdb_file])

    # get feature Matrix for Helices of PDB file
    XHelix, YHelix = fc.construct_window_matrix(helix_features,
                                                conf.helix_assigner,
                                                conf.helix_window_size)

    XStrand, YStrand = fc.construct_window_matrix(strand_features,
                                                  conf.strand_assigner,
                                                  conf.strand_window_size)

    # if we cannot extract features from this PDB file, return None
    if not XHelix:
        return None, None

    with open(pdb_file, 'r') as f:
        struc = PDBParser().get_structure(get_id(pdb_file), f)

    # predict Helix position
    helix_predictor = conf.helix_predictor
    pred_helix = helix_predictor.predict(XHelix)
    pred_helix = correct_helices(pred_helix)

    # predict Strand Positions
    strand_predictor = conf.strand_predictor
    pred_strand = strand_predictor.predict(XStrand)

    # expand predicted helix and strand annotations
    pred_helix_expand = expand(struc, pred_helix,
                               conf.helix_window_size)
    pred_strand_expand = expand(struc, pred_strand,
                                conf.strand_window_size)

    true_helix_expand = expand(struc, YHelix,
                               conf.helix_window_size)
    true_strand_expand = expand(struc, YStrand,
                                conf.strand_window_size)

    # merge separate predictions of helices and sheets
    pred_expand = merge_prediction(pred_helix_expand, pred_strand_expand)
    true_expand = merge_prediction(true_helix_expand, true_strand_expand)

    return true_expand, pred_expand
