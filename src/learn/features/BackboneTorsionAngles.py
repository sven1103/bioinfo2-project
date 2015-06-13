from Bio.PDB import PDBParser
import numpy as np

from src.learn.features.WindowFeature import WindowFeature
from src.learn.features.SheetFeature import SheetFeature
from src.pdb.extract import get_backbone_torsion_angles, get_amino_acids
from src.util import get_id


class BackboneTorsionAngles(WindowFeature, SheetFeature):

    def __init__(self):
        super(BackboneTorsionAngles, self).__init__()
        self.torsion_angles = None

    def tell_context(self, pdb_path):

        with open(pdb_path, 'r') as f:
            struc = PDBParser().get_structure(get_id(pdb_path), f)

        self.torsion_angles = get_backbone_torsion_angles(get_amino_acids(struc))

    def encode(self, entity):

        return [f(pos)
                for pos in map(lambda x: x.get_id()[1],  entity)
                for f in (lambda x: self.torsion_angles[x][0],
                          lambda x: self.torsion_angles[x][1])]

    def encode_sheet(self, strand1, strand2):

        strand1_phi = []
        strand1_psi = []
        strand2_phi = []
        strand2_psi = []

        for pos in strand1:

            angles = self.torsion_angles[pos]
            strand1_phi.append(angles[0])
            strand1_psi.append(angles[1])

        for pos in strand2:
            angles = self.torsion_angles[pos]
            strand2_phi.append(angles[0])
            strand2_psi.append(angles[1])

        return [np.mean(strand1_phi), np.mean(strand1_psi),
                np.mean(strand2_phi), np.mean(strand2_psi)]
