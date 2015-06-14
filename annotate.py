import argparse
import sys
import os

from src.conf import conf
from src.learn.base import annotate
from src.learn.target_encoding import Q3_MAPPING, TYPE_MAPPING
from src.util import get_id

PARSER_DESCRIPTION = ("Annotates PDB file with secondary structure elements "
                      "and prints final annotation to several files. "
                      "The Annotation includes Helix Positions, Strand "
                      "Positions, and Sheets together with their orientation.")

def strand_annotate(pos, sheets):

    for (sheet_id, sheet) in enumerate(sheets):
        for strand in sheet:

            if strand[0] <= pos <= strand[1]:

                return sheet_id, strand[2]

    return ' ', ' '


def main(argv):

    # configure command-line parser for the annotation main
    parser = argparse.ArgumentParser(description=PARSER_DESCRIPTION)
    parser.add_argument('PDB', type=str, help="PDB file to be annotated.")
    parser.add_argument('-o', type=str,
                        help="Where the output files should go to.",
                        default='out')

    # parse command-line arguments
    args = parser.parse_args(argv[1:])

    # check whether PDB file actually exists
    if not os.path.isfile(args.PDB):
        sys.stderr.write('No PDB file: ' + args.PDB)
        sys.exit(1)

    # create output directory, if it does not yet exist
    if not os.path.exists(args.o):
        os.makedirs(args.o)

    conf.set_dir('.')
    conf.load_predictors()

    # get PDBID
    pdbid = get_id(args.PDB)

    # annotate given PDB file
    (true_expand, pred_expand, true_sheets, pred_sheets) = annotate(args.PDB)

    with open(args.o + os.path.sep + pdbid + '_predicted', 'w') as f:

        for (pos, cls) in pred_expand.items():

            sheet_id, orient = strand_annotate(pos, pred_sheets)

            f.write('\t'.join(map(str, [pos,cls,Q3_MAPPING[cls],
                                        TYPE_MAPPING[cls], sheet_id, orient]))
                    + os.linesep)

    with open(args.o + os.path.sep + pdbid + '_true', 'w') as f:
        for (pos, cls) in true_expand.items():

            sheet_id, orient = strand_annotate(pos, true_sheets)

            f.write('\t'.join(map(str, [pos,cls,Q3_MAPPING[cls],
                                        TYPE_MAPPING[cls],sheet_id,orient]))
                    + os.linesep)

if __name__ == '__main__':
    main(sys.argv)











