# bioinfo2-project

The following scripts make up the entry points for the ML-based secondary structure annotation:

annotate.py:

Annotates a given PDB file with its secondary structure and produces two Output files:

 * PDBID_true: The secondary structure annotation as it is given by the PDB file
 * PDBID_predicted: The predicted secondary structure.

Error Codes:
1: PDB file is not present at the specified path


create_predictors.py:

Trains HELIX, STRAND, and SHEET prediction models using all PDB files as training data
specified by the argument. Make sure that you set features, you want to use, in feature_list.py
and the parameters of the predictors in configuration.py.


evaluate.py

Script that is intended to produce evaluation statistics of current trained predictors.

Error Codes:
3: Provided CV method is not available.



########################################################################################################

pdb.hydrogen.py:

Given a Bio.PDB.Structure.Structure object struc (using PDBParser.get_structure() for instance) you can
now simply add all missing backbone hydrogens with

sha = StructureHydrogenAdder(struc)
sha.supplement()
struc_hydrogen = sha.struc

Then struc_hydrogen is a Structure object containing all missing hydrogens of the backbone.














