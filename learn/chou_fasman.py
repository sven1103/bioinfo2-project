__author__ = 'lukas'

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


class ChouFasmanHelix(object):

    def __init__(self):
        pass

    def encode(self, entity):
        """
        Encodes list of residues `entity` with corresponding `Chou-Fasman`
        parameters.

        :param entity: Sequence of Residue objects which should be Amino Acids.
        :return: Chou Fasman encoding of this entity as list
        """
        return [HELICES[aa.get_resname().upper()] for aa in entity]
