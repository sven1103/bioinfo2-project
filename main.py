__author__ = 'lukas'
import argparse
import sys


PARSER_DESCRIPTION = "Predicts secondary structure annotation for given PDB file"



def main(argv):

    parser = argparse.ArgumentParser(description=PARSER_DESCRIPTION)
    parser.add_argument('PDB', type=str, help="PDB file to be annotated.")
    parser.add_argument('')

    # get Predictor for basic Secondary Structure annotation













if  __name__ == '__main__':
    main(sys.argv)
