import os

from Bio.PDB import PDBParser
import numpy as np

from pdb.hydrogen import _compute_hydrogen, _validate
from util import window, angle
import pdb.extract as extract

training_path = 'material/training/'
training = os.listdir(training_path)

nmr = 0

fout = open("angles", 'w')

for f in map(lambda x: training_path + x, training):

    struc = PDBParser().get_structure("test", f)
    # fetch amino acids of file
    for (aa1, aa2) in window(extract.get_amino_acids(struc)):

        if not (_validate(aa1) and _validate(aa2)):
            continue

        atoms2 = aa2.get_unpacked_list()

        # extract required atoms
        nitrogen = np.array(atoms2[0].get_coord())
        ca = np.array(atoms2[1].get_coord())

        # look for hydrogen in residue:
        hydrogen = _compute_hydrogen(aa1, aa2)

        # If there is no annotated backbone hydrogen, skip this residue
        if hydrogen is None:
            continue

        vec1 = hydrogen - nitrogen
        vec2 = ca - nitrogen

        fout.write(str(angle(vec1, vec2, deg=True)) + os.linesep)

fout.close()
