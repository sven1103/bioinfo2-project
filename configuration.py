from sklearn.ensemble import RandomForestClassifier


#############################################################################
# If you want to train new predictors, declare them here
##############################################################################
helix_predictor = RandomForestClassifier(n_estimators=20)
strand_predictor = RandomForestClassifier(n_estimators=20)
sheet_predicor = RandomForestClassifier(n_estimators=20)


helix_window_size = 5
strand_window_size = 3
temp_dir = 'temp'
pred_dir = 'predictors'
eval_dir = 'evaluation'
