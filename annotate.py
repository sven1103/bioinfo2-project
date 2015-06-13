import argparse
import sys
import os

from src.conf import conf

PARSER_DESCRIPTION = ("Annotates PDB file with secondary structure elements "
                      "and prints final annotation to several files. "
                      "The Annotation includes Helix Positions, Strand "
                      "Positions, and Sheets together with their orientation.")

def main(argv):

    # configure command-line parser for the annotation main
    parser = argparse.ArgumentParser(description=PARSER_DESCRIPTION)
    parser.add_argument('PDB', type=str, help="PDB file to be annotated.")

    # parse command-line arguments
    args = parser.parse_args(argv[1:])

    # load Predictors
    conf.load_predictors(os.path.dirname(__file__) + os.path.sep +
                         'predictors')






if __name__ == '__main__':
    main(sys.argv)
