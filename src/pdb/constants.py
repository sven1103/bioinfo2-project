AMINO_ACIDS = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
               "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
               "TYR", "VAL"]

NMR = 'nmr'

Q1 = 0.42
Q2 = 0.20
DIMENSIONAL_FACTOR = 332
HBOND_THRESHOLD = -0.5  # kcal/mole
DERIVATION_H = 1e-3
HNCA_ANGLE_DEG = 119
NH_DISTANCE2 = 0.9409
NH_DISTANCE = 0.97

# Useful regular expressions.

# Matches Sequence of whitespace with any length (at least one char)
RE_WHITESPACE = r'[ \t\n\r\f\v]+'

# Matches any PDBID
RE_PDBID = r'[1-9][0-9a-zA-Z]{3}'


# Matches file that contains hydrogen bond data
RE_HB = r'.*\.hb$'

# Matches PDB files
RE_PDB = r'.*\.pdb'
