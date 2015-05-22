"""
Provides utility functions that do not fit anywhere else
"""
from itertools import islice

import numpy as np


def window(seq, n=2):
    """Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ..."""
    it = iter(seq)
    result = tuple(islice(it, n))

    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def angle(vec1, vec2, deg=True):
    """
    Computes angle between vec1 and vec2. Make
    sure that the vectors have the same dimension.

    :param vec1:  First vector.
    :param vec2:  Second vector.
    :param deg: Whether to return angle in degrees instead of radians.
    :return: Angle between both vectors
    """
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)

    res = np.arccos((np.dot(vec1, vec2)) / (norm1 * norm2))

    return res if not deg else np.degrees(res)
