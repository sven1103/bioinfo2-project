from pdb.extract import get_amino_acids
from util import window


class WindowExtractor(object):

    def __init__(self,
                 struc,
                 window_size,
                 features,
                 positions):
        """
        Identify each AA segment of length `window_size` as point in
        Feature Space and extract specified features for all training points
        of `struc`

        :param struc: BioPython Structure Object
        :param window_size: Number of consequtive amino acids that should be
        considered as one entity.
        """
        self.struc = struc
        self.positions = positions
        self.features = features
        self.window_size = window_size

    def _encode(self, entity):

        # TODO There should be a better, more pythonic way, to accomplish this
        encoding = []
        for feature in self.features:
            encoding.extend(feature.encode(entity))

        return encoding

    def entities(self):

        for entity in window(get_amino_acids(self.struc, self.positions),
                             self.window_size):

            yield self._encode(entity)
