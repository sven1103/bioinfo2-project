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
