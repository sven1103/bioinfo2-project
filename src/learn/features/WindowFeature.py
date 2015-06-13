"""
All Features that want to encode Windows of consecutive amino acid sequences
must inherit from this class.
"""

class WindowFeature(object):

    def __init__(self):
        pass

    def encode(self, entity):
        pass

    def tell_context(self, pdb_path):
        return None
