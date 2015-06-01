"""
Extract some required attributes from PDB files using
the BioPython library.
"""
import math
import warnings

from Bio import BiopythonWarning, PDB

import constants


# remove the ugly warning from the PDB module, that chains are discontinuous
warnings.simplefilter("ignore", BiopythonWarning)


def get_amino_acids(struc):
    """
    Gets all amino acids of specified PDB file
    :param struc: Structure objects as returned by BioPythons PDB file parser
    :return: Generator of all Amino Acid residues of PDB file
    """
    for residue in struc.get_residues():
        if residue.get_resname() in constants.AMINO_ACIDS:
            yield residue

def get_atom_max_serial_number(struc):

    # Get list of all atoms and map them to their serial numbers,
    # then return max
    return max(map(lambda atom: atom.get_serial_number(),
                   struc.get_atoms()))


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


def compute_torsion_angles(residue, next_residue):
    """
    Little helper function, calculates the backbone phi and psi torsion
    angles from the given residues and returns them
    :param residue: The amino acid residue the torsion angles shall be computed
    :return: Phi and psi backbone torsion angles
    """
    # extract the atoms for the torsion calculation
    # 1.) for the phi
    atom_CO_1 = residue['C'].get_vector()
    atom_N_1 = residue['N'].get_vector()
    atom_CA_1 = residue['CA'].get_vector()
    atom_CA_2 = next_residue['CA'].get_vector()
    atom_CO_2 = next_residue['C'].get_vector()
    atom_N_2 = next_residue['N'].get_vector()

    phi_angle = PDB.calc_dihedral(atom_CO_1, atom_N_2, atom_CA_2, atom_CO_2)
    psi_angle = PDB.calc_dihedral(atom_N_1, atom_CA_1, atom_CO_1, atom_N_2)

    # convert into degrees
    return math.degrees(phi_angle), math.degrees(psi_angle)


def get_backbone_torsion_angles(generator_aa, pos_helix, pos_sheet):
    """
    Calculates the backbone torsion angles for annotated alpha-helices and
    beta-sheets from a Generator with all amino acid residues from a PDB-file
    :param generator_aa: Generator with extracted amino acid residues from the
    PDB-file
    :return: Two lists of tuples, containing the phi and psi torsion angles
    for alpha-helices and beta-sheets
    """
    torsion_angles_helix = []
    torsion_angles_sheet = []

    for residue in generator_aa:
        residue_pos = residue.get_id()[1]
        chain = residue.get_full_id()[2]
        if chain is not 'A':
            break
        try:
            next_residue = generator_aa.next()
        except StopIteration:
            break
        if residue_pos in pos_helix and next_residue.get_id()[1] in pos_helix:
            torsion_angles_helix.append(compute_torsion_angles(residue,
                                                               next_residue))
        elif residue_pos in pos_sheet and next_residue.get_id()[1] in pos_sheet:
            torsion_angles_sheet.append(compute_torsion_angles(residue,
                                                               next_residue))
    return torsion_angles_helix, torsion_angles_sheet


if __name__ == "__main__":
    # test file
    pdb_file = "/home/sven/Git/bioinformatics2/assignment_2/pdb/3VSU.pdb"
    #pdb_file = "/home/fillinger/git/bioinformatics2/assignment_2/pdb/3VSU.pdb"
    struct = PDB.PDBParser().get_structure("test", pdb_file)
    residues = get_amino_acids(struct)
    # example: print the residues that are in helices
    structure_positions = get_secondary_structure_annotation(pdb_file)
    print structure_positions[0]
    res = get_backbone_torsion_angles(residues, structure_positions[0],
                                      structure_positions[1])
    # alpha helices torsions
    for angles in res[0]:
        print angles[0], ":", angles[1]

    # beta sheet torsions
    print "-----------------"
    for angles in res[1]:
        print angles[0], ":", angles[1]