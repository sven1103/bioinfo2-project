__author__ = 'fillinger'

import target_encoding
from Bio import PDB
from src.pdb import extract as ex
import difflib

class PatternAnnotator(object):

    def __init__(self, hb_file, pdb_file):
        self.hb_file = open(hb_file, "r").read().splitlines()
        self.protein = PDB.PDBParser().get_structure("protein", pdb_file)
        self.pdb_file = pdb_file
        self.structure_patterns = self.find_patterns()
        self.annotation = get_pdb_annotation(self.protein, self.pdb_file)
        self.annotation_prediction = self.refine_structure_annotation()

    def find_patterns(self):
        list_annotation = []
        list_pairs = []
        for row in self.hb_file:

            pair = (int(row.split(" ")[0]), int(row.split(" ")[1]))
            list_pairs.append(pair)
            list_annotation.append(self.estimate_structure_type(pair))
        return zip(list_pairs, list_annotation)

    def estimate_structure_type(self, hb_pair):
        distance = abs(hb_pair[0] - hb_pair[1])

        if distance == 3:
            type = target_encoding.TARGET_CODES.get('310-Helix')
        elif distance == 4:
            type = target_encoding.TARGET_CODES.get('Alpha-Helix')
        elif distance == 5:
            type = target_encoding.TARGET_CODES.get('Pi-Helix')
        elif distance > 5:
            type = target_encoding.TARGET_CODES.get('sss')
        else:
            type = 0
        return type

    def refine_structure_annotation(self):

        list_structure_part = []
        helix_object = []
        strand_object = []

        for residue in self.structure_patterns:
            if str(residue[1]) in "123":
                if strand_object:
                    list_structure_part.append(strand_object)
                    strand_object = []
                helix_object.append(residue)
            elif str(residue[1]) in "5":
                if helix_object:
                    list_structure_part.append(helix_object)
                    helix_object = []
                strand_object.append(residue)

        list_structure_part.append(helix_object)
        list_structure_part.append(strand_object)

        refined_list_helices = []
        refined_list_sheets = []

        """
        Now eliminate residues from the annotation
        """
        for line in list_structure_part:
            for elem in line:
                if str(elem[1]) in "5":
                    if has_neighbor(elem, list_structure_part, elem[1]):
                        refined_list_sheets.append(elem)
                elif str(elem[1]) in "123":
                    if has_neighbor(elem, list_structure_part, elem[1]):
                        refined_list_helices.append(elem)

        helices = group_helices(refined_list_helices)
        sheets = group_sheets(refined_list_sheets)

        """
        Make the complete annotation
        """
        print sheets
        all_residues = zip(*self.annotation)[0]
        annotation_pred = []
        for residue in all_residues:
            if residue in helices:
                annotation_pred.append("H")
            elif residue in sheets:
                annotation_pred.append("S")
            else:
                annotation_pred.append("C")

        print all_residues
        print annotation_pred


        refined_annotation = remove_gaps("".join(annotation_pred), 3)

        return zip(all_residues, refined_annotation)

    def accuracy_prediction(self):
        annotation_real = "".join(zip(*self.annotation)[1])
        annotation_pred = "".join(zip(*self.annotation_prediction)[1])


        print annotation_pred
        print annotation_real


        counter_common_structure = 0
        counter = 0

        for res in annotation_real:
            if annotation_pred[counter] is res:
                counter_common_structure += 1
            counter += 1

        return float(counter_common_structure)/len(annotation_real)


def has_neighbor(res, list_annotations, type):

    if str(type) in "123":
        for line in list_annotations:
            for element in line:
                if str(element[1]) in "123":
                    if element[0][0] == res[0][0] - 1\
                            or element[0][0] == res[0][0] + 1:
                        return True
    elif str(type) in "5":
        for line in list_annotations:
            for element in line:
                if str(element[1]) in "5":
                    if element[0][0] == res[0][0] - 2 \
                            or element[0][0] == res[0][0] + 2:
                        return True
    return False


def remove_gaps(string_annotation, window_size):

    counter = 0
    window_list = []

    while counter in range(0, len(string_annotation)):
        window_string = string_annotation[counter:counter+window_size]
        if window_string is "SCS":
            window_list.append(["SSS"])
        elif window_string is "HCH":
            window_list.append(["HHH"])
        else:
            window_list.append(window_string)
        counter += 1

    print window_list
    refined_annotation = []
    for window in window_list:
        refined_annotation.append(window[0])

    print "".join(refined_annotation)

    return refined_annotation







def group_helices(list_helix_residues):
    list_helices = []
    helix_object = []

    for res in list_helix_residues:
        if not helix_object:
            helix_object.append(res)
        else:
            if res[0][0] == helix_object[len(helix_object)-1][0][0] + 1:
                helix_object.append(res)
            else:
                if len(helix_object) > 1:
                    list_helices.append(helix_object)
                helix_object = []
    list_helices.append(helix_object)

    final = set([])
    for line in list_helices:
        if len(line) > 1:
            for elem in line:
                final.add(elem[0][0])
                final.add(elem[0][1])
    return final


def group_sheets(list_sheet_residues):
    list_sheets = set([])
    sheet_object = []
    final_list = []
    for res in list_sheet_residues:
        for compare_res in list_sheet_residues:
            if compare_res[0][0] == res[0][0] + 2:
                list_sheets.add(res)
                list_sheets.add(compare_res)
                break

    for res in sorted(list_sheets):
        if not sheet_object:
            sheet_object.append(res)
        else:
            if res[0][0] == sheet_object[len(sheet_object)-1][0][0] + 2:
                sheet_object.append(res)
            else:
                final_list.append(sheet_object)
                sheet_object = []
                sheet_object.append(res)

    final = set([])
    for line in final_list:
        if len(line) > 1:
            for elem in line:
                final.add(elem[0][0])
                final.add(elem[0][1])
    return list(final)


def get_pdb_annotation(protein_structure, pdb_file):
    positions_helices, position_sheets, _ = ex.get_secondary_structure_annotation(pdb_file)

    all_aa = ex.get_amino_acids(protein_structure)

    combine_helices = []

    for key in positions_helices.keys():
        combine_helices.extend(positions_helices[key])

    annotation = []
    residues = []

    for residue in all_aa:
        residues.append(residue.get_id()[1])
        if residue.get_id()[1] in combine_helices:
            annotation.append("H")
        elif residue.get_id()[1] in position_sheets:
            annotation.append("S")
        else:
            annotation.append("C")

    return zip(residues, annotation)


if __name__ == "__main__":
    hb_file = "../training_data2/1q4k.hb"
    pdb_path = "/home/fillinger/git/bioinformatics2/assignment_2/pdb/1q4k.pdb"
    hbond_pattern = PatternAnnotator(hb_file, pdb_path)
    #hbond_pattern.refine_structure_annotation()
    hbond_pattern.accuracy_prediction()