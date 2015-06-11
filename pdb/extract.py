"""
Extract some required attributes from PDB files using
the BioPython library.
"""

__author__ = 'fillinger'
import math
import warnings
from collections import defaultdict

from Bio import BiopythonWarning, PDB

import constants






# remove the ugly warning from the PDB module that chains are discontinuous
warnings.simplefilter("ignore", BiopythonWarning)


def get_amino_acids(struc, positions=None, chain='A'):
    """
    Gets all amino acids of specified PDB file
    :param struc: Structure objects as returned by BioPythons PDB file parser
    :return: Generator of all Amino Acid residues of PDB file
    """
    for residue in struc.get_residues():

        segid = residue.get_id()[1]

        # skip if this residue does not belong to the chain of interest.
        if chain != residue.get_full_id()[2]:
            continue

        if residue.get_resname() in constants.AMINO_ACIDS:

            if positions is not None and segid not in positions:
                continue
            else:
                yield residue


def get_atom_max_serial_number(struc):

    # Get list of all atoms and map them to their serial numbers,
    # then return max
    return max(map(lambda atom: atom.get_serial_number(),
                   struc.get_atoms()))


# TODO Add possibility to extract Sheet annotation

def get_secondary_structure_annotation(path):
    """
    Reads a PDB-file and returns two lists of residue positions in normal
    right handed helices and beta sheets. We need the residue positions and its
    annotated structure for calculating torsion angles for only these amino
    acids.
    Unfortunately, the PDB-module from BioPython does not support secondary
    structure read out so far...
    :param path: the path to the pdb-file
    :return: a tuple of two lists of residue positions in [0]: helices
                                                          [1]: beta sheets
    """
    # A list with residue positions in sheets
    record_strand_aa = []
    # A dictionary with helix types and corresponding positions
    record_helix_aa = defaultdict(list)

    # keep track of start and end positions of all encountered sheets
    record_sheets = []
    current_sheet = []
    strand_counter = 0

    with open(path, 'r') as pdb_file:
        for line in pdb_file:
            # if line provides a HELIX information
            if line.startswith("HELIX") and "A" in line[19]:
                # we want only the residues of one chain to avoid duplicates
                # extract helix class

                hclass = int(line[39:40])
                # record the aminoacid residue position in a helix
                record_helix_aa[hclass].extend(range(int(line[21:25]),
                                                     int(line[33:37]) + 1))

            # if line provides a SHEET information
            elif line.startswith("SHEET") and "A" in line[21]:
                # if line provides a SHEET information
                start_res = int(line[22:26])  # the start residue of a sheet
                term_res = int(line[33:37])  # the terminating res. of a sheet
                record_strand_aa += range(start_res, term_res, 1)

                # extract strand type, assume 0 if we encounter somewhat
                # invalid entry

                try:
                    sclass = int(line[38:40])
                except ValueError:
                    sclass = 0

                current_sheet.append((start_res, term_res, sclass))

                strand_counter += 1

                # extract number of sheets in strands
                strand_number = int(line[14:16])

                # if this is the last strand of this sheet, we reset the
                # the strand  and add this sheet to all sheets
                if strand_counter == strand_number:

                    strand_counter = 0
                    record_sheets.append(list(current_sheet))
                    current_sheet = []

    return record_helix_aa, record_strand_aa, record_sheets


def compute_torsion_angles(previous_residue, residue, next_residue):
    """
    Little helper function, calculates the backbone phi and psi torsion
    angles from the given residues and returns them
    :param residue: The amino acid residue the torsion angles shall be computed
    :return: Phi and psi backbone torsion angles
    """
    # print previous_residue.get_id()[1], residue.get_id()[1], next_residue.get_id()[1]
    # extract the atoms for the torsion calculation
    # 1.) for the phi
    atom_CO_0 = previous_residue['C'].get_vector()
    atom_N_1 = residue['N'].get_vector()
    atom_CA_1 = residue['CA'].get_vector()
    atom_CO_1 = residue['C'].get_vector()
    atom_N_2 = next_residue['N'].get_vector()

    phi_angle = PDB.calc_dihedral(atom_CO_0, atom_N_1, atom_CA_1, atom_CO_1)
    psi_angle = PDB.calc_dihedral(atom_N_1, atom_CA_1, atom_CO_1, atom_N_2)

    # convert into degrees
    return math.degrees(phi_angle), math.degrees(psi_angle)


def get_backbone_torsion_angles(generator_aa, pos_helix=None, pos_sheet=None):
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

    if pos_helix is None:
        previous_residue = 0
        for residue in generator_aa:
            if previous_residue == 0:
                # print "init"
                previous_residue = residue
                curr_residue = residue
                next_residue = generator_aa.next()
            else:
                previous_residue = curr_residue
                curr_residue = next_residue
                next_residue = residue

            chain = residue.get_full_id()[2]
            if chain is not 'A':
                break
            torsion_angles_helix.append(
                compute_torsion_angles(previous_residue,
                                       curr_residue,
                                       next_residue))
    else:
        previous_residue = 0
        for residue in generator_aa:
            if previous_residue == 0:
                previous_residue = residue
                curr_residue = residue
                next_residue = generator_aa.next()
            else:
                previous_residue = curr_residue
                curr_residue = next_residue
                next_residue = residue

            chain = residue.get_full_id()[2]
            if chain is not 'A':
                break

            if curr_residue.get_id()[1] in pos_helix and\
                    next_residue.get_id()[1] in pos_helix:
                torsion_angles_helix.append(
                    compute_torsion_angles(previous_residue,
                                           curr_residue,
                                           next_residue))
            elif curr_residue.get_id()[1] in pos_sheet and\
                    next_residue.get_id()[1] in pos_sheet:
                torsion_angles_sheet.append(
                    compute_torsion_angles(previous_residue,
                                           curr_residue,
                                           next_residue))
    return torsion_angles_helix, torsion_angles_sheet


if __name__ == "__main__":
    # test file
    pdb_file = "/home/sven/Git/bioinformatics2/assignment_2/pdb/1SMC.pdb"
    #pdb_file = "/home/fillinger/git/bioinformatics2/assignment_2/pdb/3PSD.pdb"
    struct = PDB.PDBParser().get_structure("test", pdb_file)
    residues = get_amino_acids(struct)
    # example: print the residues that are in helices
    structure_positions = get_secondary_structure_annotation(pdb_file)
    print structure_positions[0][1]
    # res = get_backbone_torsion_angles(residues, structure_positions[0][1],
    #                                  structure_positions[1])
    res = get_backbone_torsion_angles(residues)
    # alpha helices torsions
    for angles in res[0]:
        print angles[0], ":", angles[1]

    # beta sheet torsions
    print "-----------------"
    for angles in res[1]:
        print angles[0], ":", angles[1]
