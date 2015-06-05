__author__ = 'fillinger'


""" The encoding of the targets during the machine learning process.
As far as i know numbers can be arbitrary chosen, as long as they are unique
and integer?
"""
TARGET_CODES = {"x-s-s": 1, "s-s-s": 2, "s-s-x": 3, "x-x-x": -1}

# class encoding used for helix prediction
HELIX_CODES = {'No Helix': 0, 'Alpha': 1, 'Pi': 2, '310': 3}
