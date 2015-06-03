__author__ = 'fillinger'


import numpy as np

class FeatureMatrix:
    """
    This class will provide a feature matrix object, containing the features
    and the targets as preparation for the machine learning process.
    This class should ensure, that feature space and targets have the same length.
    """
    def __init__(self, file_list):
        self.file_list = file_list
        self.feature_matrix, self.target = self.build_features_from_file()

    def build_features_from_file(self):
        for file in self.file_list:
            print "nothing"
        return -1,-1




