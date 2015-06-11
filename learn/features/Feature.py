"""
Just a generic abstraction from a Feature that we want to consider.
"""

class Feature(object):

    def __init__(self):
        self.name = 'Feature'

    def encode(self, entity):
        pass

    def tell_context(self, pdb_path):
        return None
