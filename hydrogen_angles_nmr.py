__author__ = 'lukas'
import os

from Bio.PDB import PDBParser
import numpy as np

import pdb.constants as co
from pdb.property import experimental_method
from util import window, angle
import pdb.extract as extract

training_path = 'material/training/'
training = os.listdir(training_path)

nmr = 0


fout = open("angles_nmr", 'w')

for f in map(lambda x: training_path + x, training):

    # only consider NMR files.
    if co.NMR in experimental_method(f):

        struc = PDBParser().get_structure("test", f)
        # fetch amino acids of file
        for (aa1, aa2) in window(extract.get_amino_acids(struc)):

            atoms2 = aa2.get_unpacked_list()

            # extract required atoms
            nitrogen = np.array(atoms2[0].get_coord())
            ca = np.array(atoms2[1].get_coord())

            # look for hydrogen in residue:
            hydrogen = None
            for atom in atoms2:
                if atom.get_name().strip() == "H":
                    hydrogen = np.array(atom.get_coord())
                    break

            # If there is no annotated backbone hydrogen, skip this residue
            if hydrogen is None:
                continue

            vec1 = ca - nitrogen
            vec2 = hydrogen - nitrogen

            fout.write(str(angle(vec1, vec2, deg=True)) + '\n')

fout.close()
