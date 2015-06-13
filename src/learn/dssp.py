from __future__ import division
import src.pdb.constants as co


def potential(r_ON, r_CH, r_OH, r_CN):

    return co.Q1 * co.Q2 * ((1/r_ON) + (1/r_CH) - (1/r_OH) - (1/r_CN)) *\
        co.DIMENSIONAL_FACTOR
