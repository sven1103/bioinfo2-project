__author__ = 'lukas'

from Bio.PDB import PDBParser

from src.pdb.extract import get_secondary_structure_annotation
from src.util import get_id, window
from src.learn.WindowExtractor import WindowExtractor
from sheets import sheet_encode
from src.pdb.hydrogen import StructureHydrogenAdder


class FeatureContext(object):
    """
    Allows the construction of feature matrices by defining features
    used for encoding entities.
    """

    def __init__(self, pdb_files):

        # extract secondary structure annotations from all PDB files, so
        # map each filepath to the secondary structures
        self.sec_struc = \
            {pdb_file: get_secondary_structure_annotation(pdb_file)
             for pdb_file in pdb_files}

    def _update(self, pdb_path, features):
        """
        Sets all features into the new context of the PDB file.

        :param pdb_path: Path to the PDB files that defines the context
        :param features: List of features as defined by the "features" package
        """
        for feature in features:
            feature.tell_context(pdb_path)

    def construct_window_matrix(self, features, annotator, window_size):
        """
        Constructs feature matrix using the features and
        wanted target annotation.

        :param features: List of features as defined by the
        feature module.
        :param annotator:
        :param window_size: Window size to be used to encode the features with.
        :return: X,Y, where X encodes an entity per row with columns
        representing the features to be used.
        """

        X = []
        Y = []

        for pdb in self.sec_struc:

            # set all features to the the current PDB contest
            self._update(pdb, features)

            # set up Window Extractor for current PDB file
            with open(pdb, 'r') as f:

                struc = PDBParser().get_structure(get_id(pdb), f)

                # ensure that the current structure contains Hydrogens
                sha = StructureHydrogenAdder(struc)
                sha.supplement()
                struc = sha.struc

                helix_aa, strand_aa, _ = self.sec_struc[pdb]
                we = WindowExtractor(struc, window_size, features)

                for (positions, training_point) in we.entities():

                    X.append(training_point)
                    Y.append(annotator(positions, helix_aa, strand_aa))

        return X, Y

    def construct_sheet_matrix(self, features):

        X = []
        Y = []

        for pdb in self.sec_struc:

            # set all features to the current PDB contest
            self._update(pdb, features)

            _, _, sheets = self.sec_struc[pdb]

            # consider each sheet
            for sheet in sheets:

                # extract all linked pairs of Strands
                for (strand1, strand2) in window(sheet, 2):

                    # use all features to encode bots strands
                    training_points = sheet_encode(range(strand1[0], strand1[1] + 1),
                                                   range(strand2[0], strand2[1] + 1),
                                                   features, pdb)
                    X.append(training_points)
                    Y.append(strand2[2])
        return X, Y

    def get_sheets(self):
        """
        Return only the sheets of this Feature Context

        :return: Map from PDB path to sheets
        """
        return {k: v[2] for k, v in self.sec_struc.iteritems()}


