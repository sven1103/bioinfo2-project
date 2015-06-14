from src.learn.features.BackboneTorsionAngles import BackboneTorsionAngles
from src.learn.features.ChouFasmanHelix import ChouFasmanHelix
from src.learn.features.HydrogenBondPattern import HydrogenBondPattern


hydrogen_bonds = HydrogenBondPattern(2)
chou_fasman_helix = ChouFasmanHelix()
backbone_torsion = BackboneTorsionAngles()

helix_features = [hydrogen_bonds]
strand_features = [backbone_torsion]

sheet_features = [hydrogen_bonds, backbone_torsion]
