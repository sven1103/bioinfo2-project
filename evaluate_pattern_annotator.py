__author__ = 'sven'


from src.learn import PatternAnnotator
import os
import re
from src.pdb.constants import RE_PDB, RE_HB


def evaluate_accuracy(pdb_dir, hb_dir):
    """
    Evaluates the accuracy for secondary structure prediction of one or more
    proteins.
    :param pdb_path: path to the pdb files (*.pdb)
    :param hb_path: path to the hydrogen bond files (*.hb)
    :return: list with accuracies in per cent
    """

    list_accuracies = []
    index = 0

    # sort files
    pdb_files = sorted(filter(lambda x: re.match(RE_PDB, x) is not None,
                              os.listdir(pdb_dir)))

    hb_files = sorted(filter(lambda x: re.match(RE_HB, x) is not None,
                             os.listdir(hb_dir)))

    hb_files_filtered = filter(lambda x: x.split(".")[0] + ".pdb" in pdb_files,
                               hb_files)

    # build paths
    pdb_paths = map(lambda x: pdb_dir + os.path.sep + x, pdb_files)
    hb_paths = map(lambda x: hb_dir + os.path.sep + x, hb_files_filtered)

    # annotate every structure and evaluate accuracy
    while index < len(pdb_files):
        print "----------------------------"
        print "Evaluate ", pdb_paths[index]
        print "Evaluate ", hb_paths[index]
        try:
            protein_annotation = PatternAnnotator.PatternAnnotator(
                hb_paths[index], pdb_paths[index]
            )
            list_accuracies.append(protein_annotation.accuracy_prediction())
        except IndexError:
            print "IndexError in file ", hb_files[index]
            pass
        print "Accuracy: ", list_accuracies[len(list_accuracies)-1]
        index += 1

    print list_accuracies


if __name__ == "__main__":
    path_pdb = "training_data/pdb_files"
    path_hb = "training_data/hb_files"
    evaluate_accuracy(path_pdb, path_hb)