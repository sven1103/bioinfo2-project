"""
This file is intended to allow the specification  of features to
use for encoding AA residues. You can use all features from
src.learn.features. You may also implement own features. Then make sure to
subclass `Window Featgure` or `Sheet Feature`.
"""
# import your features here
from src.learn.features.BackboneTorsionAngles import BackboneTorsionAngles
from src.learn.features.ChouFasmanHelix import ChouFasmanHelix,\
    ChouFasmanStrand
from src.learn.features.HydrogenBondPattern import HydrogenBondPattern


# make instances of your Features by suppling required feature parameters
hydrogen_bonds = HydrogenBondPattern(2)
chou_fasman_helix = ChouFasmanHelix()
chou_fasman_strand = ChouFasmanStrand()
backbone_torsion = BackboneTorsionAngles()


# set Helix, Strand, and Sheet Features as lists of previous defined instances.
# Make sure that helix, and strand features subclass `WindowFeature` and
# that sheet features subclass `SheetFeature`
helix_features = [backbone_torsion]
strand_features = [backbone_torsion]
sheet_features = [hydrogen_bonds, backbone_torsion]
