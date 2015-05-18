"""
Extract some required attributes from PDB files using
the BioPython library.
"""
from Bio import PDB

from . import constants


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
    Unfortunately, the PDB-modul from BioPython does not support secondary
    structure read out so far...
    :param file: the path to the pdb-file
    :return: a touple of two lists of residue positions in [0]: helices
                                                           [1]: beta sheets
    """
    return -1, -1


def get_backbone_torsion_angles(generator_aa):
    """
    Calculates the backbone torsion angles for annotated alpha-helices and
    beta-sheets from a Generator with all amino acid residues from a PDB-file
    :param generator_aa: Generator with extracted amino acid residues from the
    PDB-file
    :return: Two list of touples, containing the phi and psi torsion angles for
    alpha-helices
    """
    current_residue = generator_aa.next()
    return -1
