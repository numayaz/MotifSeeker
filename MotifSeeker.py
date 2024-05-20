import argparse
from argparse import RawTextHelpFormatter

# Test
# For testing, remember to `pip install argparser` in your terminal
# To run MotifSeeker in terminal, `python MotifSeeker.py -h`



# Create parser #
parser = argparse.ArgumentParser(description= 'Welcome to MotifSeeker!\n', formatter_class=RawTextHelpFormatter)



# Create subparsers #




# Define arguments here #
# parser.add_argument('-return', type=int, help="Test return an int")
    
# Parse arguments here #
args = parser.parse_args()
    
    