import argparse
import sys

PARSER_DESCRIPTION = ("Annotates PDB file with secondary structure elements "
                      "and prints final annotation to several files. "
                      "The Annotation includes Helix Positions, Strand "
                      "Positions, and Sheets together with their orientation.")

def main(argv):

    parser = argparse.ArgumentParser(description=PARSER_DESCRIPTION)
    parser.add_argument('PDB', type=str, help="PDB file to be annotated.")

    args = parser.parse_args(argv[1:])

























if __name__ == '__main__':
    main(sys.argv)
