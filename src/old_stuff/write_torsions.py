from src.util import reformat_list

__author__ = 'sven'

import sys
import os
import csv

from Bio import PDB

from src import util
import src.pdb.extract as ex


def write_csv_angles(path_output, torsion_angles, suffix):
    suffix += ".csv"
    with open(path_output + "torsion_angles_" + suffix, "wb") as csvfile:
        writer = csv.writer(csvfile, delimiter=",")
        for element in torsion_angles:
            # print element
            writer.writerow([element[0], element[1]])


def main(args):
    torsion_angles_helices = []
    torsion_angles_sheets = []
    path_to_pdb, path_output = util.check_arguments(args)
    file_list = os.listdir(path_to_pdb)

    for protein in file_list:
        print "Reading from ", path_to_pdb + protein
        structure = PDB.PDBParser().get_structure("", path_to_pdb + protein)
        temp_pos_helices, temp_pos_sheets =\
            ex.get_secondary_structure_annotation(path_to_pdb + protein)
        residues = ex.get_amino_acids(structure)
        temp_angles_helices, temp_angles_sheets = \
            ex.get_backbone_torsion_angles(residues)
        torsion_angles_helices.append(temp_angles_helices)
        torsion_angles_sheets.append(temp_angles_sheets)

    torsion_angles_sheets = reformat_list(torsion_angles_sheets)
    torsion_angles_helices = reformat_list(torsion_angles_helices)
    write_csv_angles(path_output, torsion_angles_sheets, "sheets_all")
    write_csv_angles(path_output, torsion_angles_helices, "helices_all")

main(sys.argv[1:])