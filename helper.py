__author__ = 'sven'

import getopt
import sys


def printCommandLineHelp():
    print "\nBioinformatics II Project - Secondary Structure Prediction\n\
Version: 0.1\n\
Authors: Lukas Zimmermann and Sven Fillinger\n\
Required Python Version: 2.7\n\
Supported format: PDB\n\n\
Syntax: python2.7 main.py [options] \n\n\
Commands:\n\n\
 -h \t\t get help :)\n\
 -p [path/]  You have to provide the path to the pdb files\n" \
          " -o [path/]  output path for the calculated energies and torsions"


def check_arguments(argv):
    try:
        opts, args = getopt.getopt(argv, "hp:o:")
    except getopt.GetoptError:
        print "\n***[Error]: Wrong number of arguments.\n"
        print "***[HINT]: Did you provide the path to the pdb files?"
        printCommandLineHelp()
        sys.exit()
    if not opts:
        printCommandLineHelp()
        sys.exit()
    for opt, arg in opts:
        if opt == "-h":
            printCommandLineHelp()
            sys.exit()
        elif opt == "-p":
            path_to_pdb = arg
        elif opt == "-o":
            path_output = arg
        elif not opt:
            print "ojee"
    return path_to_pdb, path_output

