import cPickle
import os

from sklearn.ensemble import RandomForestClassifier

from evaluation import train_test_split
import configuration as conf
from learn.FeatureContext import FeatureContext

train, test = train_test_split(conf.pdb_files)
fc = FeatureContext(train)

helix_predictor = RandomForestClassifier(n_estimators=20)
strand_predictor = RandomForestClassifier(n_estimators=20)


X_train, Y_train = fc.construct_matrix(conf.helix_features,
                                       conf.helix_assigner,
                                       conf.helix_window_size)
helix_predictor.fit(X_train, Y_train)

X_train, Y_train = fc.construct_matrix(conf.strand_features,
                                       conf.strand_assigner,
                                       conf.strand_window_size)
strand_predictor.fit(X_train, Y_train)


with open(conf.predictors_dir + os.path.sep + 'random_forest_20_helix', 'w') as f:
    cPickle.dump(helix_predictor, f)

with open(conf.predictors_dir + os.path.sep + 'random_forest_20_strand', 'w') as f:
    cPickle.dump(strand_predictor, f)
