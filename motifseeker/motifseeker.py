#!/usr/bin/env python

"""
Command-line script to perform motif enrichment analysis

Similar to HOMER's findMotifsGenome.pl
"""

import argparse
import pyfaidx
from . import myutils as myutils
from mypileup import __version__
from argparse import RawTextHelpFormatter
import os
import sys
import numpy as np
import math
import seqlogo
import scipy
import random
from bed_reader import open_bed, sample_file
import sys

# Here's how HOMER performs MEF (http://homer.ucsd.edu/homer/ngs/peakMotifs.html)

# Variables/Objects required
nucs = {"A": 0, "C": 1, "G": 2, "T": 3}


def main(): 
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

    parser.add_argument("--version", help="Print the version and quit", \
		action="version", version = '{version}'.format(version=__version__))

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

    # Load FASTA
	if args.fasta_ref is not None:
		if not os.path.exists(args.fasta_ref):
			myutils.ERROR("{fasta} does not exist".format(fasta=args.fasta_ref))
		reffasta = pyfaidx.Fasta(args.fasta_ref)
	else:
		myutils.ERROR("Fasta file not provided")

    # Load and parse BED using bed_reader
    filename = sample_file(args.inputfile)
    bed = open_bed(filename)
    val = bed.read()

    # Test
    print(val.shape)
    del bed
