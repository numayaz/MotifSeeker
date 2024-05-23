# Imports #
# argparse: Python command line interface
# numpy, bed_reader: Parse .bed file

import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import numpy as np
from bed_reader import open_bed, sample_file

# `pip install argparser` `pip install numpy` `pip install bed-reader[samples,sparse]`
# To run MotifSeeker in terminal, `python MotifSeeker.py -h`
# Here's how HOMER performs MEF (http://homer.ucsd.edu/homer/ngs/peakMotifs.html)

# Create parser
parser = argparse.ArgumentParser(prog="MotifSeeker", description= "Command-line tool to perform Motif Enrichment Analysis", formatter_class=RawTextHelpFormatter)

# Input file, we want it to work on .bed files #
parser.add_argument("inputfile", help=".bed file", type=str)

# Output
parser.add_argument("-o", "--out", help="Write output to file." \
                    "Default: stdout", metavar="FILE", type=str, required=False)

# Other arguments
parser.add_argument("-r", "--ref", help="fasta reference genome file", \
                    metavar="FILE", type=str, required=False)

# Parse arguments here
args = parser.parse_args()
    
# Setup output file
if args.out is None:
    outf = sys.stdout
else: 
    outf = open(args.out, "w")

# Load and parse BED using bed_reader
filename = sample_file(args.inputfile)
bed = open_bed(filename)
val = bed.read()

# Test
print(val.shape)
del bed
