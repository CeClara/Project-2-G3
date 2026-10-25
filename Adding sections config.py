#I have edited some more sections together ( Clara DW didnt really change the code just made the varibles consistent )
# ============================================================
# PROJECT 2: ORF FINDER TOOL
# File: config.py
# Purpose: Shared imports and project settings


###STUFF to still add 
#nested partial reading frames 
#csv
#find orf corrdinates 
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

# SECTION 1A: IMPORTS

import csv
import statistics
from collections import Counter
from pathlib import Path


# AMINO ACID CODON TABLE

codon_table = {
    'GGU':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly',
    'GCU':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
    'UUU':'Phe','UUC':'Phe','UUA':'Leu','UUG':'Leu',
    'CUU':'Leu','CUC':'Leu','CUA':'Leu','CUG':'Leu',
    'AUU':'Ile','AUA':'Ile','AUC':'Ile','AUG':'Met',
    'GUU':'Val','GUC':'Val','GUA':'Val','GUG':'Val',
    'UCU':'Ser','UCC':'Ser','UCA':'Ser','UCG':'Ser',
    'AGU':'Ser','AGC':'Ser','CCU':'Pro','CCC':'Pro',
    'CCA':'Pro','CCG':'Pro','ACU':'Thr','ACC':'Thr',
    'ACA':'Thr','ACG':'Thr','UAU':'Tyr','UAC':'Tyr',
    'CAU':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
    'AAU':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys',
    'GAU':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
    'UGU':'Cys','UGC':'Cys','UGG':'Trp','CGU':'Arg',
    'CGC':'Arg','CGA':'Arg','CGG':'Arg','AGA':'Arg',
    'AGG':'Arg','UAA':'Stop','UGA':'Stop','UAG':'Stop'
}


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

VALID_DNA_BASES = ["A", "C", "G", "T"]

# REPLACED: Temporarily use 0 to find short ORFs when testing.
# Change back to 90 when analysing longer sequences.
MIN_ORF_LENGTH = 0

DEFAULT_OUTPUT_FILE = "orf_results.csv"


# FUNCTION: validate_dna()
# Removes whitespace, converts DNA to uppercase and checks
# for invalid nucleotides. Returns cleaned DNA or False.

def validate_dna(dnaseq):

    # Remove whitespace and convert to uppercase
    clean_seq = "".join(dnaseq.split()).upper()

    # Reject empty sequences
    if not clean_seq:
        return False

    # Check that every character is a valid DNA base
    for base in clean_seq:
        if base not in VALID_DNA_BASES:
            return False

    return clean_seq


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

    # Match each base with its complementary nucleotide
    complementary_seq = "".join(ReverseCompDNA[base] for base in seq)

    # Reverse the complementary strand
    return complementary_seq[::-1]


# FUNCTION: translateseq()
# Returns DNA starting at a selected frame offset.
# This does NOT translate DNA into amino acids.
# init_pos = 0, 1 or 2.

def translateseq(seq, init_pos=0):

    return seq[init_pos:]


# FUNCTION: gen_reading_frames()
# Generates all six reading frames:
# three on the forward strand and
# three on the reverse-complement strand.

def gen_reading_frames(seq):

    frames = []

    # Generate reverse complement once
    reverse_seq = reverse_complement(seq)

    # Three forward reading frames
    frames.append(translateseq(seq, 0))
    frames.append(translateseq(seq, 1))
    frames.append(translateseq(seq, 2))

    # Three reverse reading frames
    frames.append(translateseq(reverse_seq, 0))
    frames.append(translateseq(reverse_seq, 1))
    frames.append(translateseq(reverse_seq, 2))

    return frames


# ADDED: FUNCTION translate_orf()   PITA 
# Converts a DNA ORF into amino acids using the existing
# RNA codon_table dictionary.

def translate_orf(ORFsequence):

    # ADDED: Convert DNA to RNA by replacing T with U.
    RNAsequence = ORFsequence.replace("T", "U")

    aasequence = ""

    # ADDED: Read RNA three bases at a time.
    for position in range(0, len(RNAsequence) - 2, 3):

        currentCodon = RNAsequence[position:position+3]

        # ADDED: Find the corresponding amino acid.
        amino_acid = codon_table[currentCodon]

        # ADDED: Stop codons do not code for amino acids.
        if amino_acid == "Stop":
            break

        # ADDED: Add the amino acid to the protein sequence.
        aasequence += amino_acid + "-"

    # ADDED: Remove the final hyphen.
    return aasequence.rstrip("-")


# FUNCTION: searchCodon()
# Searches one reading frame for complete ORFs.
# An ORF starts with ATG and ends at the first
# in-frame stop codon: TAA, TAG or TGA.

def searchCodon(sequence, length):

    MAXLIMIT = length

    STARTcodon = START_CODONS[0]
    STOPcodon = STOP_CODONS

    position = 0
    maxORF = []

    # Read the DNA three nucleotides at a time
    while position < len(sequence) - 2:

        currentCodon = sequence[position:position+3]

        # Check whether the codon is a start codon
        if currentCodon == STARTcodon:

            start_position = position
            scan_position = position + 3

            # Search for a stop codon in the same frame
            while scan_position < len(sequence) - 2:

                currentCodon = sequence[scan_position:scan_position+3]

                if currentCodon in STOPcodon:

                    # Extract DNA from start to stop codon
                    ORFsequence = sequence[start_position:scan_position+3]

                    # Only keep ORFs meeting minimum length
                    if len(ORFsequence) >= MAXLIMIT:
                        maxORF.append(ORFsequence)

                    # REPLACED: Removed position = scan_position.
                    # This allows the outer loop to check additional
                    # ATG codons inside an ORF (nested ORFs).
                    break

                # Move to the next codon
                scan_position += 3

        # Check the next codon
        position += 3

    return maxORF


# MAIN PROGRAM
# Ask for DNA until valid input is entered.

while True:

    dnaseq = input("Enter your DNA sequence: ")

    clean_seq = validate_dna(dnaseq)

    if clean_seq != False:
        print("Valid DNA sequence:", clean_seq)
        break

    # REPLACED: Corrected the invalid-format message.
    print("Invalid format, please only use ATGC bases.")


# Generate reading frames and find ORFs.

ORFList = []
iteration = 0

frame_names = ["+1", "+2", "+3", "-1", "-2", "-3"]

print("\nReading Frames")

for frame in gen_reading_frames(clean_seq):

    print(f"\n******** FRAME {frame_names[iteration]} ********")

    ORF = searchCodon(frame, MIN_ORF_LENGTH)

    ORFList.append(ORF)

    if ORF:
        for sequence in ORF:

            print("ORF found:", sequence)
            print("Length:", len(sequence), "bases")

            # ADDED: Translate each detected ORF.
            aasequence = translate_orf(sequence)

            # ADDED: Display the amino acid sequence.
            print("Amino acid sequence:", aasequence)

    else:
        print("No complete ORFs found.")

    iteration += 1


# Find the longest ORF across all six frames.

all_orfs = []

for frame_orfs in ORFList:
    all_orfs.extend(frame_orfs)

if all_orfs:

    longest_orf = max(all_orfs, key=len)

    print("\n******** LONGEST ORF ********")
    print("DNA sequence:", longest_orf)
    print("Length:", len(longest_orf), "bases")

    # ADDED: Display the translation of the longest ORF.    IDK hope it wqorks didnt use ai for this 
    print("Amino acid sequence:", translate_orf(longest_orf))

else:
    print("\nNo complete ORFs found.")
