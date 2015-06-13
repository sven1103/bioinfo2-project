"""
Does not really encode an entity but rather just returns
the amino acids as entity.
"""
from src.learn.features import WindowFeature


class WindowIdentity(WindowFeature):

    def __init__(self):
        super(WindowIdentity, self).__init__()

    def tell_context(self, pdb_path):
        pass

    def encode(self, entity):

        return list(entity)
