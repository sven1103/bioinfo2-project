__author__ = 'lukas'
from learn.features.Feature import Feature
from pdb.extract import get_backbone_torsion_angles

class BackboneTorsionAngles(Feature):

    def __init__(self):
        super(BackboneTorsionAngles, self).__init__()

    def set_context(self, pdb_path):
        pass

    def encode(self, entity):
        res = []
        all_angles, _ = get_backbone_torsion_angles(iter(entity))

        for (phi, psi) in all_angles:
            res.append(phi)
            res.append(psi)

        return res
