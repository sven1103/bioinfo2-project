__author__ = 'lukas'

from Bio.PDB import PDBParser

from pdb.extract import get_secondary_structure_annotation
from util import get_id
from learn.WindowExtractor import WindowExtractor


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

    def construct_matrix(self, features, annotator, window_size):
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

                helix_aa, strand_aa, sheet_aa = self.sec_struc[pdb]
                we = WindowExtractor(struc, window_size, features)

                for (positions, training_point) in we.entities():

                    X.append(training_point)
                    Y.append(annotator(positions, helix_aa, strand_aa))

        return X, Y

    def get_sheets(self):
        """
        Return only the sheets of this Feature Context

        :return: Map from PDB path to sheets
        """
        return {k: v[2] for k, v in self.sec_struc.iteritems()}
