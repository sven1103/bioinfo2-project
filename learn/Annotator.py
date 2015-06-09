__author__ = 'lukas'
from target_encoding import TARGET_CODES



class Annotator(object):

    def __init__(self, classes):

        self.classes = classes

    def __call__(self, positions, helix_aa, sheet_aa):

        y = 0

        # try to annotate a beta strand
        if len(positions) == 4:

            triplet_type = ''

            for pos in positions[1:4]:

                if pos in sheet_aa:
                    triplet_type += 's'
                else:
                    triplet_type += 'x'

            if triplet_type in TARGET_CODES:
                y = TARGET_CODES[triplet_type]

        # try to annotate helix type
        helix_aa_alpha = set(helix_aa[1])
        helix_aa_310 = set(helix_aa[5])
        positions_set = set(positions)

        if positions_set <= helix_aa_alpha:
            y = TARGET_CODES['Alpha-Helix']

        elif positions_set <= helix_aa_310:
            y = TARGET_CODES['310-Helix']

        # check whether we are actually interested in this class
        if y in self.classes:
            return y

        # if not, just say we found a coil
        return TARGET_CODES['Coil']
