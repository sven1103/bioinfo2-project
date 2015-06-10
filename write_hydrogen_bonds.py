"""
Creates file for each pdb file in the data set denoting its Hydrogen Bonds.
"""
import os

from Bio.PDB import PDBParser

from pdb.hydrogen import StructureHydrogenAdder
from pdb.extract import get_amino_acids
from learn.features import HydrogenBondPattern

# directory where the hydrogen bonds should be written to
hb_dir = 'training_data2'

# where the PDB files are that should be used
# pdb_files_path = 'material/training/'
pdb_files_path = '/home/fillinger/git/bioinformatics2/assignment_2/pdb/'

# You do not need to modify anything below here


# list of all files
pdb_paths = map(lambda x: pdb_files_path + x, os.listdir(pdb_files_path))

# create a directory for hydrogen bonds, if it does not yet exist
if not os.path.exists(hb_dir):
    os.makedirs(hb_dir)

# open each file for reading
for pdb_path in pdb_paths:

    pdb_name = pdb_path.split('.')[0].split('/')[7]

    with open(pdb_path, 'r') as f:

        # open file where we want to write the hydrogen bonds into
        with open(hb_dir + os.path.sep + pdb_name + '.hb', 'w') as fw:

            struc = PDBParser().get_structure(pdb_name, f)

            # Add missing hydrogen
            sha = StructureHydrogenAdder(struc)
            sha.supplement()

            # define that HBBonds should be extracted_
            hb_pattern = HydrogenBondPattern.HydrogenBondPattern(2)

            # extract all hydrogen bonds
            for (aa1, aa2, pot1, pot2) in hb_pattern.encode(get_amino_acids(sha.struc)):

                line = str(aa1.get_id()[1]) + ' ' + str(aa2.get_id()[1]) + ' ' + str(pot1) + ' ' + str(pot2)
                fw.write(line + os.linesep)
