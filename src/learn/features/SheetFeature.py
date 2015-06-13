"""
Just a generic abstraction from a Feature that we want to consider.
"""

class SheetFeature(object):

    def __init__(self):
        pass

    def encode_sheet(self, strand1, strand2):
        pass

    def tell_context(self, pdb_path):
        return None
