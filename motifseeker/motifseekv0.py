import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import numpy as np
from . import functions

nucs = {"A": 0, "C": 1, "G": 2, "T": 3}

# Create parser
parser = argparse.ArgumentParser(prog="MotifSeeker", description= "Command-line tool to perform Motif Enrichment Analysis", formatter_class=RawTextHelpFormatter)

# Input file, we want it to work on .bed files #
parser.add_argument("inputfile", help=".bed file containing peaks", type=str)

parser.add_argument("-f", "--fasta-ref", help="fasta reference genome file", \
                        metavar="FILE", type=str)

# Other arguments
parser.add_argument("-o", "--out", help="Write output to file." \
                        "Default: stdout", metavar="FILE", type=str, required=False)
    
parser.add_argument("-p", "--pval", help="p value threshold for motif sequences" \
                         "Default: 0.05", metavar="NUMERICAL P VALUE", type=float, required=False)
    
parser.add_argument("-s", "--size", help="size of motif" \
                        "Default: 100", metavar="NUMERICAL SIZE VALUE", type=int, required=False)

# Parse arguments here
args = parser.parse_args()

if args.size is not None:
        size = args.size
else:
        size = 100
        
# Setup output file
if args.out is None:
        outf = sys.stdout
else: 
        outf = open(args.out, "w")
