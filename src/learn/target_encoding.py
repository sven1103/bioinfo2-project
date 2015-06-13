__author__ = 'fillinger'


""" The encoding of the targets during the machine learning process.
As far as i know numbers can be arbitrary chosen, as long as they are unique
and integer?
"""
TARGET_CODES = {'Coil': 0,
                'Alpha-Helix': 1,
                'Pi-Helix': 2,
                '310-Helix': 3,
                'xss': 4,
                'sss': 5,
                'ssx': 6}

Q3_MAPPING = {0: '-', 1: 'H', 2: 'H', 3: 'H',
              4: 'E', 5: 'E', 6: 'E'}

TYPE_MAPPING = {0: 'Coil', 1: 'ALPHA', 2: 'PI', 3: '310',
                4: 'STRAND', 5: 'STRAND', 6: 'STRAND'}
