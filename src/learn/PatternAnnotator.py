__author__ = 'fillinger'

import target_encoding
from Bio import PDB
from src.pdb import extract as ex
import difflib

class PatternAnnotator(object):
    """
    This class expects a pdb-file and a file containing its computed
    hydrogen bonds (below a set energy threshold, i.e. -0.5kcal/mol).
    It extracts the secondary structure annotation from the pdb-file and
    tries to predict an annotation from the hydrogen bond file. Both
    annotations can then be compared and the accuracy computed.

    In principle the annotation prediction has 5 steps so far:

        1.) Estimation of secondary structure type by hydrogen bond residue
        distance

        2.) Grouping of residues in helices and sheets

        3.) Forming concatenated structure objects within the helices and
        sheets group (search for logical sequence patterns)

        4.) Make the preliminary annotation

        5.) Refine the annotation by removing gaps in annotated structures
    """

    def __init__(self, hb_file, pdb_file):
        # Read the hydrogen bond files (*.hb)
        self.hb_file = open(hb_file, "r").read().splitlines()

        # Get the protein structure from the pdb-file
        self.protein = PDB.PDBParser().get_structure("protein", pdb_file)

        # safe the file name
        self.pdb_file = pdb_file

        # find and assign structure types (preliminary estimation)
        self.structure_patterns = self.find_patterns()

        # get the annotation from the pdb-file
        self.annotation = get_pdb_annotation(self.protein, self.pdb_file)

        # make own structure prediction by refining the annotation
        self.annotation_prediction = self.refine_structure_annotation()

    def find_patterns(self):
        """
        Extract all the hbond pairs and preliminary estimate the structure type
        """
        list_annotation = []
        list_pairs = []
        for row in self.hb_file:

            pair = (int(row.split(" ")[0]), int(row.split(" ")[1]))
            list_pairs.append(pair)
            list_annotation.append(self.estimate_structure_type(pair))
        return zip(list_pairs, list_annotation)

    def estimate_structure_type(self, hb_pair):
        """
        Receives an hydrogen bond pair and estimates the structure type
        from the residue distance.
        For i+3 we predict 310-helices, for i+4 alpha-helices and for
        i+5 pi-helices (rare). For patterns found with residue distances
        greater 5 we suggest beta sheet structures.
        """
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
        """
        From the preliminary structure estimation we calculate groups
        of concatenated helix and sheet structure elements. This will be
        further refined then.
        """
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
        Now eliminate residues from the annotation.
        We don't want single annotations without neighboring structure
        elements.
        """
        for line in list_structure_part:
            for elem in line:
                if str(elem[1]) in "5":
                    if has_neighbor(elem, list_structure_part, elem[1]):
                        refined_list_sheets.append(elem)
                elif str(elem[1]) in "123":
                    if has_neighbor(elem, list_structure_part, elem[1]):
                        refined_list_helices.append(elem)

        """
        Group the residues if they are consecutive.
        """
        helices = group_helices(refined_list_helices)
        sheets = group_sheets(refined_list_sheets)

        """
        Make the complete annotation
        """
        all_residues = zip(*self.annotation)[0]
        annotation_pred = []
        for residue in all_residues:
            if residue in helices:
                annotation_pred.append("H")
            elif residue in sheets:
                annotation_pred.append("S")
            else:
                annotation_pred.append("C")

        """
        Remove gaps in secondary structure annotation (it is unlikely,
        that a helix has one single residue in itself that is not part of a
        helix, same for sheets.
        """
        refined_annotation = remove_gaps("".join(annotation_pred), 3)

        return zip(all_residues, refined_annotation)

    def accuracy_prediction(self):
        """
        Calculates the accuracy for the predicted structures.
        :return:
        """
        annotation_real = "".join(zip(*self.annotation)[1])
        annotation_pred = "".join(zip(*self.annotation_prediction)[1])

        print "Annotation PDB:\t\t", annotation_real
        print "Annotation pred:\t", annotation_pred
        counter_common_structure = 0
        counter = 0
        counter_common_sheet = 0
        counter_common_helix = 0

        for res in annotation_real:
            if annotation_pred[counter] is res:
                counter_common_structure += 1
                if res is "H" and annotation_pred[counter] is "H":
                    counter_common_helix += 1
                if res is "S" and annotation_pred[counter] is "S":
                    counter_common_sheet += 1
            counter += 1

        acc_q3 = float(counter_common_structure)/len(annotation_real)
        try:
            acc_helix = float(counter_common_helix)/annotation_real.count('H')
        except ZeroDivisionError:
            acc_helix = "NA"
            pass
        try:
            acc_sheet = float(counter_common_sheet)/annotation_real.count('S')
        except ZeroDivisionError:
            acc_sheet = "NA"
            pass

        return [acc_q3, acc_helix, acc_sheet]


def has_neighbor(res, list_annotations, type):
    """
    Checks if residue has neighboring residues that are part of the same
    structure type.
    :param res:
    :param list_annotations:
    :param type:
    :return:
    """

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
    """
    Remove gaps by sliding over the annotation and refine the sequence.
    :param string_annotation:
    :param window_size:
    :return:
    """

    counter = 0
    window_list_old = []
    window_list_refined = []

    while counter in range(0, len(string_annotation)):

        window_string = string_annotation[counter:counter+window_size]

        if window_string in "SCS":
            window_list_old.append([window_string])
            window_list_refined.append(["SSS"])
        elif window_string in "HCH":
            window_list_old.append([window_string])
            window_list_refined.append(["HHH"])
        else:
            window_list_old.append([window_string])
            window_list_refined.append([window_string])
        counter += 1

    refinement = []
    for window in window_list_refined:
        if not refinement:
            refinement.append(window[0][0:2])
        else:
            refinement.append(window[0][1])
    refinement = refinement[0:len(refinement)-1]
    return refinement


def group_helices(list_helix_residues):
    """
    Groups helices residues together, if they are consecutive.
    :param list_helix_residues:
    :return:
    """
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
    """
    Groups residues in sheets together, if they are consecutive
    :param list_sheet_residues:
    :return:
    """

    list_sheets = set([])
    sheet_object = []
    final_list = []
    for res in list_sheet_residues:
        for compare_res in list_sheet_residues:
            if compare_res[0][0] == res[0][0] + 2:
                list_sheets.add(res)
                list_sheets.add(compare_res)
                break
    # print "----"
    found_prefix = False
    # print "hhhhhh"
    for res in sorted(list_sheets):
        if not sheet_object:
            sheet_object.append(res)
        else:
            if res[0][0] == sheet_object[len(sheet_object)-1][0][0] + 2:
                sheet_object.append(res)
            else:
                found_prefix = False
                for sheet_objects in final_list:
                    if sheet_objects[len(sheet_objects)-1][0][0] + 2 == res[0][0]:
                        if sheet_objects[len(sheet_objects)-1][0][1] + 2 == res[0][1] or\
                                sheet_objects[len(sheet_objects)-1][0][1] - 2 == res[0][1]:
                            final_list[final_list.index(sheet_objects)].append(res)
                            found_prefix = True
                if not found_prefix:
                    final_list.append(sheet_object)
                    sheet_object = []
                    sheet_object.append(res)

    final_list.append(sheet_object)

    final = set([])
    for line in final_list:
        if len(line) > 1:
            for elem in line:
                final.add(elem[0][0])
                final.add(elem[0][1])
    return list(final)


def get_pdb_annotation(protein_structure, pdb_file):
    """
    Extract the annotation from the pdb-file
    :param protein_structure:
    :param pdb_file:
    :return:
    """
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
    """
    Test the shit
    """
    hb_file = "../training_data/hb_files/3F8C.hb"
    pdb_path = "../training_data/pdb_files/3F8C.pdb"
    hbond_pattern = PatternAnnotator(hb_file, pdb_path)
    #hbond_pattern.refine_structure_annotation()
    print hbond_pattern.accuracy_prediction()