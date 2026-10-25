#I have edited some more sections together ( Clara DW didnt really change the code just made the varibles consistent ) UNDER THE CURRENT CODE 
# ============================================================
# PROJECT 2: ORF FINDER TOOL
# File: config.py
# Purpose: Shared imports and project settings
# ============================================================

# SECTION 1A: IMPORTS
import csv
import statistics
from collections import Counter
from pathlib import Path


# SECTION 1B: START AND STOP CODONS

START_CODONS = ["ATG"]
STOP_CODONS = ["TAA", "TAG", "TGA"]


# SECTION 1C: ORF DATA FORMAT

# Each ORF will be stored as a dictionary with these fields:
#
# orf_id       - unique identifier for the ORF
# start        - starting coordinate in the DNA sequence
# end          - ending coordinate in the DNA sequence
# strand       - DNA strand: "+" or "-"
# frame        - reading frame: 1, 2 or 3
# length       - length of the ORF in bases, including the stop codon
# dna_sequence - DNA sequence of the ORF
#
# Coordinates are 1-based and inclusive.
#
# Example ORF:
#
# {
#     "orf_id": "ORF_1",
#     "start": 10,
#     "end": 99,
#     "strand": "+",
#     "frame": 1,
#     "length": 90,
#     "dna_sequence": "ATG...TAA"
# }


# SECTION 1D: OTHER SHARED SETTINGS

# DNA bases accepted by the input validator.
VALID_DNA_BASES = ["A", "C", "G", "T"]

# Default minimum ORF length in bases.
MIN_ORF_LENGTH = 90

# Default output filename.
DEFAULT_OUTPUT_FILE = "orf_results.csv"












''FIXED code ( please check before deleting above)'' 

# SECTION 1A: IMPORTS

import csv
import statistics
from collections import Counter
from pathlib import Path


# SECTION 1B: START AND STOP CODONS

START_CODONS = ["ATG"]
STOP_CODONS = ["TAA", "TAG", "TGA"]


# SECTION 1C: ORF DATA FORMAT

# Each ORF will eventually be stored as a dictionary containing:
#
# orf_id       - unique identifier for the ORF
# start        - starting coordinate in the DNA sequence
# end          - ending coordinate in the DNA sequence
# strand       - DNA strand: "+" or "-"
# frame        - reading frame: 1, 2 or 3
# length       - length of the ORF in bases
# dna_sequence - DNA sequence of the ORF
#
# Coordinates will be 1-based and inclusive.


# SECTION 1D: OTHER SHARED SETTINGS

# The four valid DNA bases.
VALID_DNA_BASES = ["A", "C", "G", "T"]

# Default minimum ORF length in bases.
MIN_ORF_LENGTH = 90

# Default output filename.
DEFAULT_OUTPUT_FILE = "orf_results.csv"


# FUNCTION: validate_dna()
# Removes whitespace, converts DNA to uppercase and checks for
# invalid nucleotides. Returns the cleaned DNA or False.



def validate_dna(dnaseq):

    # Remove whitespace and convert to uppercase
    clean_seq = "".join(dnaseq.split()).upper()

    # Reject empty sequences
    if not clean_seq:
        return False

    # Check that every character is a valid DNA base
    for base in clean_seq:
        if base not in ["A", "T", "G", "C"]:
            return False

    return clean_seq


# Keep asking until valid DNA is entered
while True:

    dnaseq = input("Enter your DNA sequence: ")

    clean_seq = validate_dna(dnaseq)

    if clean_seq != False:
        print("Valid DNA sequence:", clean_seq)
        break

    print("Invalid format, please only use ATGC bases.")




    # FUNCTION: load_dna()
# Opens a text file, reads its contents and validates the DNA.
# Returns the cleaned DNA or False.

def load_dna(file_path):

    with open(file_path, "r", encoding="utf-8") as file:
        dnaseq = file.read()

    return validate_dna(dnaseq)


# FUNCTION: reverse_complement()
# Generates the complementary DNA strand and reverses it,
# so the result is in the 5' to 3' direction.

ReverseCompDNA = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

def reverse_complement(seq):

    # Match each base with its complementary nucleotide.
    complementary_seq = "".join(ReverseCompDNA[base] for base in seq)

    # Reverse the complementary strand.
    return complementary_seq[::-1]


# FUNCTION: translateseq()
# CURRENT PURPOSE: Returns DNA starting at a selected frame offset.
# This does NOT translate DNA into amino acids.
# init_pos = 0, 1 or 2.

def translateseq(seq, init_pos=0):

    return seq[init_pos:]
