import cPickle
import os
import sys
import copy

from src.util import pdb_map
from src.learn.ClassAssigner import ClassAssigner
from src.learn.target_encoding import TARGET_CODES
from configuration import helix_predictor, strand_predictor, sheet_predicor,\
    helix_window_size, strand_window_size, temp_dir, pred_dir, eval_dir


class Configuration(object):

    def __init__(self):

        self.helix_window_size = helix_window_size
        self.strand_window_size = strand_window_size

        self.path = None
        self.temp_dir = None
        self.pred_dir = None
        self.eval_dir = None
        self.hb_dict = None

        self.helix_assigner = ClassAssigner([TARGET_CODES['Coil'],
                                             TARGET_CODES['Alpha-Helix'],
                                             TARGET_CODES['310-Helix']])
        self.strand_assigner = ClassAssigner([TARGET_CODES['Coil'],
                                              TARGET_CODES['sss']])
        self.helix_predictor = None
        self.strand_predictor = None
        self.sheet_predictor = None

    def load_predictors(self):

        if not os.path.exists(self.pred_dir + os.path.sep + 'HELIX'):
            sys.stderr.write('No HELIX predictor found in ' +
                             self.pred_dir + os.path.sep)
            sys.exit(2)

        # load Helix Predictor
        with open(self.pred_dir + os.path.sep + 'HELIX', 'r') as f:
            self.helix_predictor = cPickle.load(f)

        if not os.path.exists(self.pred_dir + os.path.sep + 'STRAND'):
            sys.stderr.write('No STRAND predictor found in ' +
                             self.pred_dir + os.path.sep)
            sys.exit(2)

        # load Strand Predictor
        with open(self.pred_dir + os.path.sep + 'STRAND', 'r') as f:
            self.strand_predictor = cPickle.load(f)

        if not os.path.exists(self.pred_dir + os.path.sep + 'SHEET'):
            sys.stderr.write('No SHEET predictor found in ' +
                             self.pred_dir + os.path.sep)
            sys.exit(2)

        # load Sheet orientation Predictor
        with open(self.pred_dir + os.path.sep + 'SHEET', 'r') as f:
            self.sheet_predictor = cPickle.load(f)

    def reset_predictors(self):
        self.helix_predictor = copy.copy(helix_predictor)
        self.strand_predictor = copy.copy(strand_predictor)
        self.sheet_predictor = copy.copy(sheet_predicor)

    def set_dir(self, path):
        self.path = path
        self.temp_dir = path + os.path.sep + temp_dir
        self.pred_dir = path + os.path.sep + pred_dir
        self.eval_dir = path + os.path.sep + eval_dir
        self.hb_dict = pdb_map(self.temp_dir)

