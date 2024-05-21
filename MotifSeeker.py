import argparse
from argparse import RawTextHelpFormatter
import os
import sys

# For testing, remember to `pip install argparser` in your terminal
# To run MotifSeeker in terminal, `python MotifSeeker.py -h`
# Here's how HOMER performs MEF (http://homer.ucsd.edu/homer/ngs/peakMotifs.html)

# Create parser #
parser = argparse.ArgumentParser(description= 'Welcome to MotifSeeker!\n', formatter_class=RawTextHelpFormatter)

# Input file, we want it to work on .bed files #
parser.add_argument("input", help=".bed file", type=str)


# Other arguments #
parser.add_argument("genome", help="Genome", type=str)



# Output #
parser.add_argument("-o", "--out", help="Write output to file." \
    "Default: stdout", metavar="FILE", type=str, required=False)


# Parse arguments here #
args = parser.parse_args()
    
# Setup output file
if args.out is None:
    outf = sys.stdout
else: outf = open(args.out, "w")

# Load and parse BED
with open(args.input) as f:
    for line in f:
        L = line.strip().split()
    # Do something here
