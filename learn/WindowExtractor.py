from pdb.extract import get_amino_acids
from util import window

class WindowExtractor(object):

    def __init__(self,
                 struc,
                 window_size,
                 features,
                 positions=None):
        """
        Identify each AA segment of length `window_size` as point in
        Feature Space and extract specified features for all training points
        of `struc`

        :param struc: BioPython Structure Object
        :param window_size: Number of consecutive amino acids that should be
        considered as one entity.
        """
        self.struc = struc
        self.positions = positions
        self.features = features
        self.window_size = window_size

    def _encode(self, entity):
        """
        Encodes all features with specified features
        :param entity: Sequence of consecutive AAs
        :return: Encoding of the entity
        """

        # TODO There should be a better, more pythonic way, to accomplish this
        encoding = []
        for feature in self.features:
            encoding.extend(feature.encode(entity))

        return encoding

    def entities(self):
        """
        Encoded each entity of this PDB file.

        :return: Generator of encoded entities, together with the
        positions of the residues in the first tuple-component
        """
        for entity in window(get_amino_acids(self.struc, self.positions),
                             self.window_size):

            positions = map(lambda x: x.get_id()[1],  entity)

            # all AAs of a window must be consecutive
            if not is_consecutive(positions):
                continue

            yield (positions, self._encode(entity))


def is_consecutive(positions):
    """
    Check whether positions is a range(i,j) object.

    :param positions: List of integers
    :return: Whether positions consists of consecutive positions
    """
    return positions == list(range(positions[0], positions[-1] + 1))
