"""
Extract some required attributes from PDB files using
the BioPython library.
"""
from Bio import PDB

import constants


def get_amino_acids(ident, file):
    """
    Gets all amino acids of specified PDB file
    :param ident: Identifier as String used to identify the structure
    :param file: Path or File handle to PDB file.
    :return: Generator of all Amino Acid residues of PDB file
    """

    struct = PDB.PDBParser().get_structure(ident, file)

    for res in struct.get_residues():

        if res.get_resname() in constants.AMINO_ACIDS:
            yield res


def get_secondary_structure_annotation(file):
    """
    Reads a PDB-file and returns two lists of residue positions in normal
    right handed helices and beta sheets. We need the residue positions and its
    annotated structure for calculating torsion angles for only these amino
    acids.
    Unfortunately, the PDB-module from BioPython does not support secondary
    structure read out so far...
    :param file: the path to the pdb-file
    :return: a tuple of two lists of residue positions in [0]: helices
                                                          [1]: beta sheets
    """
    pdb_file = open(file, "r")
    record_sheet_aa = []  # A list with residue positions in sheets
    record_helix_aa = []  # A list with residue positions in helices
    for line in pdb_file.readlines():
        # if line provides a HELIX information
        if "HELIX" in line[0:5] and "A" in line[19]:
            # we want only the residues of one chain to avoid duplicates
            # check for helix class
            if int(line[39:40]) == 1:
                # '1' is default alpha helix
                # record the aminoacid residue position in a helix
                record_helix_aa += range(int(line[21:25]),
                                         int(line[33:37]) + 1)
        elif "SHEET" in line[0:5] and "A" in line[21]:
            # if line provides a SHEET information
            start_res = int(line[22:26])  # the start residue of a sheet
            term_res = int(line[33:37])  # the terminating res. of a sheet
            record_sheet_aa += range(start_res, term_res, 1)
    pdb_file.close()
    return record_helix_aa, record_sheet_aa


def get_backbone_torsion_angles(generator_aa):
    """
    Calculates the backbone torsion angles for annotated alpha-helices and
    beta-sheets from a Generator with all amino acid residues from a PDB-file
    :param generator_aa: Generator with extracted amino acid residues from the
    PDB-file
    :return: Two lists of touples, containing the phi and psi torsion angles
    for alpha-helices and beta-sheets
    """
    current_residue = generator_aa.next()
    return -1


if __name__ == "__main__":
    # test file
    pdb_file = "/home/fillinger/git/bioinformatics2/assignment_2/pdb/3VSU.pdb"
    print get_secondary_structure_annotation(pdb_file)
