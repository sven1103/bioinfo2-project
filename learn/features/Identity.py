"""
Does not really encode an entity but rather just returns
the amino acids as entity.

"""
from learn.features.Feature import Feature


class Identity(Feature):

    def __init__(self):
        super(Identity, self).__init__()
        self.name = 'Identity'

    def tell_context(self, pdb_path):
        pass

    def encode(self, entity):

        return list(entity)
