import cPickle
import os

from sklearn.ensemble import RandomForestClassifier

from learn.ClassAssigner import ClassAssigner
from learn.features.BackboneTorsionAngles import BackboneTorsionAngles
from learn.features.ChouFasmanHelix import ChouFasmanHelix
from learn.features.HydrogenBondPattern import HydrogenBondPattern
from learn.target_encoding import TARGET_CODES
import variables as vars

def reset_predictors():

    global helix_predictor
    global strand_predictor

    helix_predictor = RandomForestClassifier(n_estimators=20)
    strand_predictor = RandomForestClassifier(n_estimators=20)


def get_predictor(path):
    """
    Loads Predictor object that is stored somewhere as Pickle object

    :param path:
    :return:
    """
    with open(path, 'r') as f:
        return cPickle.load(f)


hydrogen_bonds = HydrogenBondPattern(2)
chou_fasman_helix = ChouFasmanHelix()
backbone_torsion = BackboneTorsionAngles()

helix_features = [backbone_torsion]
strand_features = [backbone_torsion]

helix_assigner = ClassAssigner([TARGET_CODES['Coil'],
                               TARGET_CODES['Alpha-Helix'],
                               TARGET_CODES['310-Helix']])
strand_assigner = ClassAssigner([TARGET_CODES['Coil'],
                                TARGET_CODES['sss']])


helix_predictor = get_predictor(vars.predictors_dir + os.path.sep + 'random_forest_20_helix')
strand_predictor = get_predictor(vars.predictors_dir + os.path.sep + 'random_forest_20_strand')
