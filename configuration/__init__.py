import os
import re

from util import pdb_map
from pdb.constants import RE_HB, RE_PDB
from variables import pdb_dir, hb_dir
from machinelearning import *
from variables import *

pdb_files = map(lambda y: pdb_dir + os.path.sep + y,
                filter(lambda x: re.match(RE_PDB, x) is not None,
                       os.listdir(pdb_dir)))
hb_files = map(lambda y: hb_dir + os.path.sep + y,
               filter(lambda x: re.match(RE_HB, x) is not None,
                      os.listdir(hb_dir)))
pdb_dict = pdb_map(pdb_files)
hb_dict = pdb_map(hb_files)