
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