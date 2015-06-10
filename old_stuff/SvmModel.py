#TODO I guess we do not need this class anymore



__author__ = 'fillinger'

from sklearn import cross_validation
from sklearn import svm


class SvmModel:

    def __init__(self,
                 feature,
                 target,
                 kernel,
                 cross_validation):
        self.feature = feature
        self.target = target
        self.cross_validation = cross_validation
        self.model = svm.SVC(kernel=kernel, C=1)
        self.scores = self.compute_cv()

    def compute_cv(self):
        return cross_validation.cross_val_score(
            self.model, self.feature, self.target, cv=self.cross_validation
        )

    def print_scores(self):
        print "Scores from ", self.cross_validation, "fold cross-validation:\n"
        print self.scores


