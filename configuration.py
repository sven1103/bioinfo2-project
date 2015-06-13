from sklearn.ensemble import RandomForestClassifier

from src.learn.features.HydrogenBondPattern import HydrogenBondPattern
from src.learn.features.ChouFasmanHelix import ChouFasmanHelix
from src.learn.features.BackboneTorsionAngles import BackboneTorsionAngles

hydrogen_bonds = HydrogenBondPattern(2)
chou_fasman_helix = ChouFasmanHelix()
backbone_torsion = BackboneTorsionAngles()

helix_features = [backbone_torsion]
strand_features = [backbone_torsion]
sheet_features = [hydrogen_bonds, backbone_torsion]


helix_window_size = 5
strand_window_size = 3
temp_dir = 'temp'



#############################################################################
##### If you want to train new predictors, declare them here
##############################################################################
helix_predictor = RandomForestClassifier(n_estimators=20)
strand_predictor = RandomForestClassifier(n_estimators=20)
sheet_predicor = RandomForestClassifier(n_estimators=20)
