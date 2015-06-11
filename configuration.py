import re
import os

from pdb.constants import RE_PDB, RE_HB
from learn.features.HydrogenBondPatternFile import get_hydrogen_bond_pattern_file
from learn.features.BackboneTorsionAngles import BackboneTorsionAngles
from learn.features.ChouFasmanHelix import ChouFasmanHelix
from learn.ClassAssigner import ClassAssigner
from learn.target_encoding import TARGET_CODES



# import some fancy sklearn predictors
from sklearn.ensemble import RandomForestClassifier


##############################################################

# on which PDB files you want to train the estimators
pdb_dir = '/home/lukas/Dropbox/BI2_project/material/training/'

# additional files with hydrogen bonds to avoid recompution
hb_dir = '/home/lukas/Dropbox/BI2_project/material/features/hydrogen_bonds'


##############################################################

# set the used window sizes for helices and strands
helix_window_size = 5
strand_window_size = 3


###############################################################
# constructs files


pdb_files = map(lambda y: pdb_dir + os.path.sep + y,
                filter(lambda x: re.match(RE_PDB, x) is not None,
                       os.listdir(pdb_dir)))


hb_files = map(lambda y: hb_dir + os.path.sep + y,
               filter(lambda x: re.match(RE_HB, x) is not None,
                      os.listdir(hb_dir)))

###############################################################
# Configure here which features you want to use for prediction


hydrogen_bonds = get_hydrogen_bond_pattern_file(hb_files)()
chou_fasman_helix = ChouFasmanHelix()
backbone_torsion = BackboneTorsionAngles()


# Set Features for Helices
helix_features = [hydrogen_bonds]

# Set Features for Strands
strand_features = [backbone_torsion]

# Set Features for Sheets
# TODO


###########################################################
# Set some required class assigner


helix_assigner = ClassAssigner([TARGET_CODES['Coil'],
                               TARGET_CODES['Alpha-Helix'],
                               TARGET_CODES['310-Helix']])

strand_assigner = ClassAssigner([TARGET_CODES['Coil'],
                                TARGET_CODES['sss'],
                                TARGET_CODES['ssx'],
                                TARGET_CODES['xss']])

############################################################
# Configure Predictors for helices and sheets

helix_predictor = RandomForestClassifier(n_estimators=20)
strand_predictor = RandomForestClassifier(n_estimators=20)

