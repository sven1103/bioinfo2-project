"""
This is the central configuration file. You can set important subdirectories
and meta parameters here.
"""
from sklearn.ensemble import RandomForestClassifier

#############################################################################
# If you want to train new predictors, declare them here.
# You can use all predictors from scikit-learn.
##############################################################################
helix_predictor = RandomForestClassifier(n_estimators=20)
strand_predictor = RandomForestClassifier(n_estimators=20)
sheet_predicor = RandomForestClassifier(n_estimators=20)

# The Window size to be used for Helices
helix_window_size = 5

# The Window size to be used for Strands
strand_window_size = 3

# Directory for caching features
temp_dir = 'temp'

# Where the trained predictors should be pickled to
pred_dir = 'predictors'

# Where evaluation results should be stored
eval_dir = 'evaluation'
