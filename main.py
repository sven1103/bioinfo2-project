__author__ = 'sven'

import sys
import os
import csv

from Bio import PDB

import helper
import pdb.extract as ex


def reformat_list(list):
    """
    Reformats a list with lists of tuples into a single list of tuples
    :param list: a list of lists with touples
    :return: a single list of tuples
    """
    new_list = []
    for element in list:
        if element:
            for tuples in element:
                new_list.append(tuples)
    return new_list


def write_csv_angles(path_output, torsion_angles, suffix):
    suffix += ".csv"
    with open(path_output + "torsion_angles_" + suffix, "wb") as csvfile:
        writer = csv.writer(csvfile, delimiter=",")
        for element in torsion_angles:
            print element
            writer.writerow([element[0], element[1]])


def main(args):
    torsion_angles_helices = []
    torsion_angles_sheets = []
    path_to_pdb, path_output = helper.check_arguments(args)
    file_list = os.listdir(path_to_pdb)

    for protein in file_list:
        print "Reading from ", path_to_pdb + protein
        structure = PDB.PDBParser().get_structure("", path_to_pdb + protein)
        temp_pos_helices, temp_pos_sheets =\
            ex.get_secondary_structure_annotation(path_to_pdb + protein)
        residues = ex.get_amino_acids(structure)
        temp_angles_helices, temp_angles_sheets = \
            ex.get_backbone_torsion_angles(residues, temp_pos_helices[1],
                                           temp_pos_sheets)
        torsion_angles_helices.append(temp_angles_helices)
        torsion_angles_sheets.append(temp_angles_sheets)

    torsion_angles_sheets = reformat_list(torsion_angles_sheets)
    torsion_angles_helices = reformat_list(torsion_angles_helices)
    write_csv_angles(path_output, torsion_angles_sheets, "sheets")
    write_csv_angles(path_output, torsion_angles_helices, "helices")

main(sys.argv[1:])