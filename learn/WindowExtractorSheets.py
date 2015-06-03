__author__ = 'sven'

import pdb.extract as ex
from Bio import BiopythonWarning, PDB


class WindowExtractorSheets(object):
    """ This class will create a lists containing windows of size '3' (amino
    acids) and their respective torsion angles, as well as their classification
    as sheet or not sheet (helix or coil remaining).
    It will also distinguish between introducing sheet windows and terminating
    sheet windows and assign own classes to them.
    """
    def __init__(self,
                 struc,
                 window_size,
                 pdb_file):
        """
        WindowExtractorSheets object.
        Out of a Structure object, we will extract all information for sheets
        and non sheets respectively. We will use a window based assay and slide
        along the primary sequence and define:
            - beginning strand triple [x-s-s] -> 1
            - in strand triple [s-s-s] -> 2
            - terminating strand triple [s-s-x] -> 3
            - non strand triple [x-x-x] -> -1

        :param struc: BioPython Structure Object
        :param window_size: Number of consecutive amino acids that should be
        considered as one entity.
        :param pdb_file: The pdb file to read out
        """
        self.struc = struc
        self.window_size = window_size
        self.pdb = pdb_file

    def extract_features(self):
        """
        Extract the features from the given pdb file
        :return: a list containing the features and a list containing the
        target encoding as integer values
        """
        all_aa = list(ex.get_amino_acids(self.struc))
        _, strand_aa = ex.get_secondary_structure_annotation(self.pdb)



        print all_aa
        print strand_aa
        return [-1]


def calculate_torsions_triplet(aa_list):
    """
    Will compute the three torsion angle pairs for the respective triple
    :param aa_list:
    :return: a list with 3 tuples of torsion angles
    """
    return []


def evaluate_triplet_type(aa_list_all, aa_sheets):
    """
    Evaluates a triple and assign its type
    :param aa_list_all:
    :param aa_sheets:
    :return:
    """


if __name__ == "__main__":
    """ You can test the class and its functions here """
    pdb_test_file = "/home/sven/Git/bioinformatics2/assignment_2/pdb/3PSD.pdb"
    # make a structure object
    struc = PDB.PDBParser().get_structure("test", pdb_test_file)
    # test the class
    test = WindowExtractorSheets(struc, 3, pdb_test_file)





