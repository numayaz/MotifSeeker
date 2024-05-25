# MotifSeeker

This was created as a final project for Professor Melissa Gymrek's CSE 185 class at UCSD. MotifSeeker is a python command-line tool which aims to identify enriched motifs in the genome when provided reference genome and peak files. It aims to emulate some limited functionality of the HOMER python package.

## Installation

`MotifSeeker` requires the following python libraries to be installed:
- numpy
- argparse
- bed_reader
- scipy
- seqlogo

They can be installed using `pip`:

```pip install argparse numpy bed-reader scipy seqlogo```

### Install instructions

To install our tool, you can run the following commands:

```
git clone https://github.com/jhii7/MotifSeeker.git
cd motifseeker
python setup.py install
```

Note: if you do not have root access, you can run the commands above with additional options to install locally:

```
pip install --user argparse numpy bed-reader scipy seqlogo
python setup.py install --user
```

If the install was successful, typing `motifseeker --help` should show a useful message.

### Usage instructions

`MotifSeeker` uses the following input files from the command-line arguments:
- `.bed file` containing the peak data
- `.fa file` containing the reference genome

### Basic usage

The basic usage of `motifseeker` is:
```
motifseeker <.bed file> <.fa ref genome file> [other options]
```

To run mypileup on a small test example (using files in this repo):
```
motifseeker example_files/test.bed example_files/test.fa
```

### Additional command line options

### Input file format requirements

The BED file should have a minimum of 6 tab separated columns (additional columns will be ignored)
(Read the [HOMER documentation](http://homer.ucsd.edu/homer/ngs/peakMotifs.html#:~:text=The%20findMotifsGenome.pl%20program%20is,the%20enrichment%20of%20known%20motifs.) for more on acceptable input files)
- Column1: chromosome
- Column2: starting position
- Column3: ending position
- Column4: Unique Peak ID
- Column5: not used
- Column6: Strand (+/- or 0/1, where 0="+", 1="-")

The fasta file should have the following format
```
>chr[name]
[chromosome sequence]
```

## Contributors
