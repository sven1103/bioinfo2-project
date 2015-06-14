__author__ = 'lukas'
from src.learn.features.WindowFeature import WindowFeature

HELICES = {'GLU': 1.53,
           'ALA': 1.45,
           'LEU': 1.34,
           'HIS': 1.24,
           'MET': 1.20,
           'GLN': 1.17,
           'TRP': 1.14,
           'VAL': 1.14,
           'PHE': 1.12,
           'LYS': 1.07,
           'ILE': 1.00,
           'ASP': 0.98,
           'THR': 0.82,
           'SER': 0.79,
           'ARG': 0.79,
           'CYS': 0.77,
           'ASN': 0.73,
           'TYR': 0.61,
           'PRO': 0.59,
           'GLY': 0.53}

STRANDS = {'MET': 1.67,
           'VAL': 1.65,
           'ILE': 1.60,
           'CYS': 1.30,
           'TYR': 1.29,
           'PHE': 1.28,
           'GLN': 1.23,
           'LEU': 1.22,
           'THR': 1.20,
           'TRP': 1.19,
           'ALA': 0.93,
           'ARG': 0.90,
           'GLY': 0.81,
           'ASP': 0.80,
           'LYS': 0.74,
           'SER': 0.72,
           'HIS': 0.71,
           'ASN': 0.65,
           'PRO': 0.62,
           'GLU': 0.26}

class ChouFasmanHelix(WindowFeature):

    def __init__(self):
        super(ChouFasmanHelix, self).__init__()

    def encode(self, entity):
        """
        Encodes list of residues `entity` with corresponding `Chou-Fasman`
        parameters.

        :param entity: Sequence of Residue objects which should be Amino Acids.
        :return: Chou Fasman encoding of this entity as list
        """
        return [HELICES[aa.get_resname().upper()] for aa in entity]

class ChouFasmanStrand(WindowFeature):

    def __init__(self):
        super(ChouFasmanStrand, self).__init__()

    def encode(self, entity):
        """
        Encodes list of residues `entity` with corresponding `Chou-Fasman`
        parameters.

        :param entity: Sequence of Residue objects which should be Amino Acids.
        :return: Chou Fasman encoding of this entity as list
        """
        return [STRANDS[aa.get_resname().upper()] for aa in entity]
