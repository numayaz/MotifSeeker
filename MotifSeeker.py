import argparse
from argparse import RawTextHelpFormatter
import os


# For testing, remember to `pip install argparser` in your terminal
# To run MotifSeeker in terminal, `python MotifSeeker.py -h`



# Create parser #
parser = argparse.ArgumentParser(description= 'Welcome to MotifSeeker!\n', formatter_class=RawTextHelpFormatter)

# Input file, we want it to work on peaks.txt files #
# TODO: maybe add .bed files as input too?
parser.add_argument("peak", help="Peaks.txt file", type=str)

# Other arguments #
parser.add_argument("genome", help="Genome", type=str)



# Output #
# parser.add_argument("-o", "--out", help="Write output to file." \
#   "Default: stdout", metavar="FILE", type=str, required=False)


# Parse arguments here #
args = parser.parse_args()
    
    