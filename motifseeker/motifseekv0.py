import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import numpy as np
import pandas as pd
import scipy
from Bio import SeqIO

nucs = {"A": 0, "C": 1, "G": 2, "T": 3}

##### Define functions here: #####

# Using bed file, splice genome file and output to new FASTA file.
def READ(bedfile, genomefile):
    

    # Open bedfile
    bed = pd.read_csv(bedfile, header = None, delimiter = '\t')
    print(bed)
    beddict = bed.to_dict()
    print(beddict)

    fasta_sequences = SeqIO.parse(open(genomefile),'fasta')
"""
    with open(output_file) as out_file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            new_sequence = some_function(sequence)
            write_fasta(out_file)
"""


# Using spliced genome file, create different reading frames of varying lengths
# ex. parse7 is an array which contains reads of lengths
def PARSEGENOME(splicedGenome):
    parse7 = []
    parse8 = []
    parse9 = []
    parse10 = []
    parse11 = []
    parse12 = []
    parse13 = []
    parse14 = []
    parse15 = []
    parse16 = []
    parse17 = []
    parse18 = []
    parse19 = []
    parse20 = []
    parse21 = []
    parse22 = []
    parse23 = []

    
##################################

# Create parser
parser = argparse.ArgumentParser(prog="MotifSeeker", description= "Command-line tool to perform Motif Enrichment Analysis", formatter_class=RawTextHelpFormatter)

# Input file, we want it to work on .bed files #
parser.add_argument("inputfile", help=".bed file containing peaks", metavar ="FILE", type=str)

parser.add_argument("genome", help="fasta reference genome file", \
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

if (args.inputfile is not None):
       print("Input file loaded...")
if (args.genome is not None):
       print("Genome file loaded...")
READ(args.inputfile, args.genome)
