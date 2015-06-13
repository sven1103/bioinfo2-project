import cPickle
import os
import copy

from src.learn.ClassAssigner import ClassAssigner
from src.learn.target_encoding import TARGET_CODES
import configuration as conf


class Configuration(object):

    def __init__(self):

        self.helix_features = conf.helix_features
        self.stand_features = conf.strand_features
        self.sheet_features = conf.sheet_features

        self.helix_window_size = conf.helix_window_size
        self.strand_window_size = conf.strand_window_size

        self.temp_dir = conf.temp_dir

        self.helix_assigner = ClassAssigner([TARGET_CODES['Coil'],
                                             TARGET_CODES['Alpha-Helix'],
                                             TARGET_CODES['310-Helix']])
        self.strand_assigner = ClassAssigner([TARGET_CODES['Coil'],
                                              TARGET_CODES['sss']])
        self.helix_predictor = None
        self.strand_predictor = None
        self.sheet_predictor = None

    def load_predictors(self, path):

        # load Helix Predictor
        with open(path + os.path.sep + 'HELIX', 'r') as f:
            self.helix_predictor = cPickle.load(f)

        # load Strand Predictor
        with open(path + os.path.sep + 'STRAND', 'r') as f:
            self.strand_predictor = cPickle.load(f)

        # load Sheet orientation Predictor
        with open(path + os.path.sep + 'SHEET', 'r') as f:
            self.sheet_predictor = cPickle.load(f)

    def reset_predictors(self):
        self.helix_predictor = copy.copy(conf.helix_predictor)
        self.strand_predictor = copy.copy(conf.strand_predictor)
        self.sheet_predictor = copy.copy(conf.sheet_predicor)
