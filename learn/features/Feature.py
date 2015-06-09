"""
Just a generic abstraction from a Feature that we want to consider.
"""

class Feature(object):

    def __init__(self):
        pass

    def encode(self, entity):
        pass

    def set_context(self, pdb_path):
        pass
