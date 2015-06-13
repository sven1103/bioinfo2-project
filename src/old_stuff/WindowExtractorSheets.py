__author__ = 'sven'

from Bio import PDB

import src.pdb.extract as ex


class WindowExtractorSheets(object):
    """ This class will create a lists containing windows of size '3' (amino
    acids) and their respective torsion angles, as well as their classification
    as sheet or not sheet (helix or coil remaining).
    It will also distinguish between introducing sheet windows and terminating
    sheet windows and assign own classes to them.
    """
    def __init__(self,
                 struc,
                 pdb_file):
        """
        WindowExtractorSheets object.
        Out of a Structure object, we will extract all information for sheets
        and non sheets respectively. We will use a window based assay and slide
        along the primary sequence and define:
            - beginning strand triple [x-s-s] -> 1
            - in strand triple [s-s-s] -> 2
            - terminating strand triple [s-s-x] -> 3
            - non strand triple [x-x-x] -> 0

        :param struc: BioPython Structure Object
        :param pdb_file: The pdb file to read out
        """
        self.struc = struc
        self.pdb = pdb_file
        self.list_all_aa = []
        self.feature_list = []
        self.strand_aa = []

    def compute_features(self):
        """
        Extract the features from the given pdb file
        :return: a list containing the features and a list containing the
        target encoding as integer values
        """
        # extract all amino acids from the PDB Structure
        all_aa = ex.get_amino_acids(self.struc)
        # convert to list
        self.list_all_aa = list(all_aa)
        # get the sheet annotation
        _, self.strand_aa = ex.get_secondary_structure_annotation(self.pdb)
        # compute all torsion angles
        all_torsion_angles, _ = ex.get_backbone_torsion_angles(
            ex.get_amino_acids(self.struc)
        )

        # the window slider of size 3, slides over the primary sequence
        # and assigns position types
        for residue_pos in range(len(self.list_all_aa)):
            # the break condition, when the window is at the end
            if residue_pos+2 == len(self.list_all_aa)-1:
                break
            # build the torsion angles and the type
            self.feature_list.append((
                all_torsion_angles[residue_pos:residue_pos+3],
                self.evaluate_triplet_type(
                    self.list_all_aa[residue_pos].get_id()[1]))
            )

    def get_features(self):
        """
        Get the feature list
        :return: the feature list
        """
        return self.feature_list

    def evaluate_triplet_type(self, residue_pos):
        """
        Evaluates a triplet and assign its type
        :return:
        """
        # the string that will store the triplet type
        triplet_type = ""
        for curr_res in range(residue_pos, residue_pos+3):
            if curr_res in self.strand_aa:
                triplet_type += "s"
            else:
                triplet_type += "x"
        # assign the type
        if triplet_type in "xss":
            type_classifier = 1
        elif triplet_type in "sss":
            type_classifier = 2
        elif triplet_type in "ssx":
            type_classifier = 3
        else:
            type_classifier = 0
        return type_classifier


if __name__ == "__main__":
    """ You can test the class and its functions here """
    pdb_test_file = "/home/sven/Git/bioinformatics2/assignment_2/pdb/3PSD.pdb"
    # make a structure object
    struc = PDB.PDBParser().get_structure("test", pdb_test_file)
    # test the class
    obj = WindowExtractorSheets(struc, pdb_test_file)
    obj.compute_features()
    for torsion_set in obj.get_features():
        print torsion_set




