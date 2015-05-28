# bioinfo2-project


pdb.hydrogen.py:

Given a Bio.PDB.Structure.Structure object struc (using PDBParser.get_structure() for instance) you can
now simply add all missing backbone hydrogens with

sha = StructureHydrogenAdder(struc)
sha.supplement()
struc_hydrogen = sha.struc

Then struc_hydrogen is a Structure object containing all missing hydrogens of the backbone.














