AMINO_ACIDS = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
               "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
               "TYR", "VAL"]

NMR = 'nmr'

#UNIT_ELECTRON_CHARGE = 1.60217657e-19  # Coulomb
Q1 = 0.42  #  * UNIT_ELECTRON_CHARGE
Q2 = 0.20  #  * UNIT_ELECTRON_CHARGE
DIMENSIONAL_FACTOR = 332
HBOND_THRESHOLD = -0.5  # kcal/mole
DERIVATION_H = 1e-3
HNCA_ANGLE_DEG = 119
NH_DISTANCE2 = 0.9409
NH_DISTANCE = 0.97

# Useful regular expressions.

# Matches Sequence of whitespace with any length (at least one char)
RE_WHITESPACE = r'[ \t\n\r\f\v]+'
