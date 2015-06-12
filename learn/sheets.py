from itertools import combinations

import numpy as np

__author__ = 'lukas'


def predict_sheets(hb_pairs, strand_positions):

    # all determined strands
    strands = []

    for (a, b) in hb_pairs:

        # if the distance is not sufficient, skip this bond
        if np.abs(a - b) < 5:
            continue

        # skip this hydrogen bond, if it does not connect sheets
        if a not in strand_positions and b not in strand_positions:
            continue

        left_a = a - 1
        right_a = a + 1
        left_b = b - 1
        right_b = b + 1

        # try to extend the sheets as most as possible
        while True:

            strand1 = range(left_a, right_a + 1)
            strand2 = range(left_b, right_b + 1)

            if left_a - 1 in strand_positions:
                left_a -= 1
                continue

            if right_a + 1 in strand_positions:
                right_a += 1
                continue

            if left_b - 1 in strand_positions:
                left_b -= 1
                continue

            if right_b + 1 in strand_positions:
                right_b += 1
                continue
            break

        strands.append(strand1)
        strands.append(strand2)

    # all strands have been collected, try to merge them where possible
    merged = merge_strands(strands)
    sheets = join_strands(merged, hb_pairs)

    for sheet in sheets:
        for i in range(len(sheet)):

            a = min(sheet[i][0])
            b = max(sheet[i][0])
            sheet[i] = (a, b, sheet[i][1])

    return sheets


def connected(positions1, positions2, pairs):

    for (left, right) in [(c, d)
                        for c in positions1
                        for d in positions2]:

        if (left, right) in pairs or (right, left) in pairs:

            return True

    return False

def merge_strands(strands):

    # first, map all strands to the set interface
    strands = map(set, strands)
    merge_complete = False

    while not merge_complete:

        merge_complete = True

        for (i, j) in combinations(range(len(strands)), 2):

            if strands[i] & strands[j]\
                    or min(strands[i]) - 1 == max(strands[j])\
                    or min(strands[j]) - 1 == max(strands[i]):

                merge_complete = False

                strands.append(strands[i] | strands[j])
                del strands[j]
                del strands[i]
                break
    return strands


def join_strands(strands, hb_pairs):

    strands = sorted(strands, key=lambda x:  np.mean(list(x)))

    res = []
    strand_added = False

    for strand in strands:

        for sheet in res:

            elem = sheet[-1]
            if connected(elem[0], strand, hb_pairs):

                # TODO Predict sheet orientation here
                sheet.append((strand, 3))
                strand_added = True
                break

        if strand_added:
            strand_added = False
            continue
        # strand could not be chained to another sheet, thus we introduce
        # a new sheet here
        else:
            res.append([(strand, 0)])
    return res





