from sklearn.ensemble import RandomForestClassifier

from learn.ClassAssigner import ClassAssigner
from learn.features.BackboneTorsionAngles import BackboneTorsionAngles
from learn.features.ChouFasmanHelix import ChouFasmanHelix
from learn.features.HydrogenBondPattern import HydrogenBondPattern
from learn.target_encoding import TARGET_CODES

__author__ = 'lukas'
hydrogen_bonds = HydrogenBondPattern(2)
chou_fasman_helix = ChouFasmanHelix()
backbone_torsion = BackboneTorsionAngles()
helix_features = [hydrogen_bonds]
strand_features = [backbone_torsion]
helix_assigner = ClassAssigner([TARGET_CODES['Coil'],
                               TARGET_CODES['Alpha-Helix'],
                               TARGET_CODES['310-Helix']])
strand_assigner = ClassAssigner([TARGET_CODES['Coil'],
                                TARGET_CODES['sss'],
                                TARGET_CODES['ssx'],
                                TARGET_CODES['xss']])
helix_predictor = RandomForestClassifier(n_estimators=20)
strand_predictor = RandomForestClassifier(n_estimators=20)