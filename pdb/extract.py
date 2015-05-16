"""
Extract some required attributes from PDB files using
the BioPython library.
"""
from Bio import PDB

from . import constants


def get_amino_acids(ident, file):

    struct = PDB.PDBParser().get_structure(ident, file)

    for res in struct.get_residues():

        if res.get_resname() in constants.AMINO_ACIDS:
            yield res
