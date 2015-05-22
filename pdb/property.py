from Bio.PDB import PDBParser


def experimental_method(pdb_path):
    """
    Get String representation of Experimental method used file of interest.
    Use header for this information.

    :param pdb_path: Path to PDB file
    :return:
    """
    parser = PDBParser(get_header=True)
    parser.get_structure('', pdb_path)

    return parser.get_header()['structure_method']
