#I have edited some more sections together ( Clara DW didnt really change the code just made the varibles consistent )
# basic stuff added but some stats or some nice formatting would be nice now GUI needed to be done
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

# Each ORF is stored as a dictionary containing:
#
# orf_id       - unique identifier for the ORF
# start        - starting coordinate in the original DNA
# end          - ending coordinate in the original DNA
# strand       - DNA strand: "+" or "-"
# frame        - reading frame: 1, 2 or 3
# length       - length of the ORF in bases
# dna_sequence - DNA sequence of the ORF
#
# Original DNA coordinates will be 1-based and inclusive.


# SECTION 1D: OTHER SHARED SETTINGS

VALID_DNA_BASES = ["A", "C", "G", "T"]

# CHANGED: Length is in DNA bases, including the stop codon when present.
# Enter 0 at the prompt to disable the minimum length filter.
MIN_ORF_LENGTH = 90

DEFAULT_OUTPUT_FILE = "orf_results.csv"


# FUNCTION: validate_dna()
# Removes whitespace, converts DNA to uppercase and checks
# for invalid nucleotides. Returns cleaned DNA or False.

def validate_dna(dnaseq):

    # Remove whitespace and convert to uppercase
    clean_seq = "".join(dnaseq.split()).upper()

    # Reject empty sequences
    if not clean_seq:
        # CHANGED: Explain why empty input is invalid.
        print("Input is empty. Please enter DNA bases.")
        return False

    # Check that every character is a valid DNA base
    for base in clean_seq:
        if base not in VALID_DNA_BASES:
            # CHANGED: N is deliberately rejected rather than treated as a DNA base.
            print(f"Invalid character found: {base}. Please only use ATGC bases.")
            return False

    return clean_seq


# FUNCTION: load_dna()
# Opens a text file, reads its contents and validates the DNA.
# Returns the cleaned DNA or False.

def load_dna(file_path):

    # CHANGED: Give a friendly message if the text file cannot be opened.
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            dnaseq = file.read()
    except OSError as error:
        print(f"Could not open DNA file: {error}")
        return False

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


# FUNCTION: translate_orf()
# Converts a DNA ORF into amino acids using the existing
# RNA codon_table dictionary.

def translate_orf(ORFsequence):

    # Convert DNA to RNA by replacing T with U
    RNAsequence = ORFsequence.replace("T", "U")

    aasequence = ""

    # Read RNA three bases at a time
    for position in range(0, len(RNAsequence) - 2, 3):

        currentCodon = RNAsequence[position:position+3]

        # Find the corresponding amino acid
        amino_acid = codon_table[currentCodon]

        # Stop codons do not code for amino acids
        if amino_acid == "Stop":
            break

        # Add the amino acid to the protein sequence
        aasequence += amino_acid + "-"

    # Remove the final hyphen
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

                        # NEW: Store the ORF with its frame-relative
                        # positions instead of only its DNA sequence.
                        # These positions are 0-based, end exclusive.
                        maxORF.append({
                            "dna_sequence": ORFsequence,
                            "start_position": start_position,
                            "end_position": scan_position + 3
                        })

                    # Continue checking other ATG codons,
                    # including those inside this ORF.
                    break

                # Move to the next codon
                scan_position += 3

        # Check the next codon
        position += 3

    # CHANGED: Only complete ORFs in this SAME frame are tested for nesting.
    for ORF in maxORF:

        ORF["nested"] = False

        # NEW: Compare with every other detected ORF.
        for other_ORF in maxORF:

            if (
                other_ORF["start_position"] < ORF["start_position"]
                and other_ORF["end_position"] >= ORF["end_position"]
            ):

                # NEW: Mark this ORF as nested.
                ORF["nested"] = True
                break

    # NEW: Identify these as complete ORFs (with an in-frame start and stop).
    for ORF in maxORF:
        ORF["type"] = "complete"

    return maxORF


# NEW: Detect simplified partial-ORF candidates in one shifted reading frame.
# Inspired by partial5 / partial3 handling in urmi-21/orfipy.
# These are possible coding regions, not confirmed genes.
# A stop-delimited segment without ATG is a missing-start candidate.
# A terminal segment containing ATG is missing-stop; otherwise missing-both.
def find_partial_orfs(sequence, length):
    partialORF = []
    segment_start = 0
    end_position = len(sequence) - len(sequence) % 3

    for position in range(0, end_position, 3):
        if sequence[position:position+3] in STOP_CODONS:
            if position > segment_start and not any(
                sequence[i:i+3] in START_CODONS
                for i in range(segment_start, position, 3)
            ):
                partialORF.append({
                    "dna_sequence": sequence[segment_start:position+3],
                    "start_position": segment_start,
                    "end_position": position+3,
                    "type": "missing_start", "nested": False
                })
            segment_start = position + 3

    starts = [i for i in range(segment_start, end_position, 3)
              if sequence[i:i+3] in START_CODONS]
    for i in starts:
        partialORF.append({
            "dna_sequence": sequence[i:end_position],
            "start_position": i, "end_position": end_position,
            "type": "missing_stop", "nested": False
        })
    if not starts and end_position > segment_start:
        partialORF.append({
            "dna_sequence": sequence[segment_start:end_position],
            "start_position": segment_start, "end_position": end_position,
            "type": "missing_both", "nested": False
        })
    # NEW: Apply the same length threshold to partial candidates.
    return [ORF for ORF in partialORF if len(ORF["dna_sequence"]) >= length]


# NEW: Convert zero-based shifted-frame positions to original DNA
# 1-based inclusive coordinates. Reverse-strand start <= end.
# Inspired by coordinate mapping in 001pa/ORF-Finder.
def original_coordinates(ORF, frame_number, genome_length):
    offset = frame_number % 3
    left = ORF["start_position"] + offset
    right = ORF["end_position"] + offset
    if frame_number < 3:
        return left + 1, right
    return genome_length - right + 1, genome_length - left


# NEW: Export flat ORF dictionaries, not a nested list of six lists.
# Uses Python's built-in csv module.
def export_csv(ORFList, filename):
    columns = ["orf_id", "type", "nested", "strand", "frame",
               "start", "end", "length", "dna_sequence", "amino_acid_sequence"]
    with open(filename, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ORFList)


# MAIN PROGRAM
# Ask for DNA until valid input is entered.

while True:

    # CHANGED: Accept pasted DNA or a .txt file path.
    dnaseq = input("Enter DNA sequence or path to a .txt file: ").strip()
    if dnaseq.lower().endswith(".txt"):
        clean_seq = load_dna(dnaseq)
    else:
        clean_seq = validate_dna(dnaseq)

    if clean_seq != False:
        print("Valid DNA sequence:", clean_seq)
        break

    # CHANGED: The validator has already explained the invalid input.


# NEW: Let the user override the default length for short sequences.
while True:
    answer = input(f"Minimum ORF length in DNA bases [{MIN_ORF_LENGTH}]: ").strip()
    if not answer:
        break
    if answer.isdigit() and int(answer) >= 0:
        MIN_ORF_LENGTH = int(answer)
        # CHANGED: Zero explicitly disables the minimum-length filter.
        if MIN_ORF_LENGTH == 0:
            print("Minimum length filter disabled; very short candidates may appear.")
        break
    print("Please enter a non-negative whole number.")

# Generate reading frames and find ORFs.

ORFList = []
iteration = 0

frame_names = ["+1", "+2", "+3", "-1", "-2", "-3"]

print("\nReading Frames")

for frame in gen_reading_frames(clean_seq):

    print(f"\n******** FRAME {frame_names[iteration]} ********")

    ORF = searchCodon(frame, MIN_ORF_LENGTH)

    # CHANGED: Add candidate partial ORFs and assess nesting in this frame.
    ORF.extend(find_partial_orfs(frame, MIN_ORF_LENGTH))
    for ORF_result in ORF:
        ORF_result["nested"] = any(
            other["start_position"] <= ORF_result["start_position"]
            and other["end_position"] >= ORF_result["end_position"]
            and (other["end_position"] - other["start_position"])
            > (ORF_result["end_position"] - ORF_result["start_position"])
            for other in ORF
        )

    # NEW: Add original DNA coordinates and structured fields for CSV export.
    for ORF_result in ORF:
        ORF_result["start"], ORF_result["end"] = original_coordinates(
            ORF_result, iteration, len(clean_seq))
        ORF_result["strand"] = "+" if iteration < 3 else "-"
        ORF_result["frame"] = iteration % 3 + 1
        ORF_result["length"] = len(ORF_result["dna_sequence"])
        ORF_result["amino_acid_sequence"] = translate_orf(ORF_result["dna_sequence"])

    ORFList.append(ORF)

    if ORF:

        # NEW: Each ORF is now a dictionary containing
        # the DNA sequence, positions and nested status.
        for ORF_result in ORF:

            sequence = ORF_result["dna_sequence"]

            print("ORF found:", sequence)
            print("Length:", len(sequence), "bases")

            # CHANGED: Partial results are candidates; nesting includes partials.
            print("ORF type:", ORF_result["type"], "(candidate)" if ORF_result["type"] != "complete" else "")
            print("Nested within same frame:", ORF_result["nested"])
            # NEW: Display coordinates and strand on the original DNA.
            print("Coordinates:", ORF_result["start"], "-", ORF_result["end"])
            print("Strand:", ORF_result["strand"], "Frame:", ORF_result["frame"])

            # Translate each detected ORF
            aasequence = translate_orf(sequence)

            # Display the amino acid sequence
            # CHANGED: Partial translations are only candidate proteins.
            print("Candidate amino acid sequence:" if ORF_result["type"] != "complete" else "Amino acid sequence:", aasequence)

    else:
        print("No ORF candidates found at this length.")

    iteration += 1


# Find the longest ORF across all six frames.

all_orfs = []

for frame_orfs in ORFList:
    all_orfs.extend(frame_orfs)

# CHANGED: A partial ORF cannot be the longest COMPLETE ORF.
complete_orfs = [ORF for ORF in all_orfs if ORF["type"] == "complete"]
if complete_orfs:

    # NEW: The ORF results are dictionaries now,
    # so compare their DNA sequence lengths.
    longest_orf = max(
        complete_orfs,
        key=lambda ORF: len(ORF["dna_sequence"])
    )

    # NEW: Extract the DNA sequence from the dictionary.
    longest_sequence = longest_orf["dna_sequence"]

    print("\n******** LONGEST ORF ********")
    print("DNA sequence:", longest_sequence)
    print("Length:", len(longest_sequence), "bases")

    # NEW: Translate the longest ORF sequence.
    print("Amino acid sequence:", translate_orf(longest_sequence))

else:
    print("\nNo complete ORFs found.")

# NEW: Assign unique IDs and save CSV after all six frames are processed.
for index, ORF in enumerate(all_orfs, start=1):
    ORF["orf_id"] = f"ORF_{index}"

# NEW: This is a console program; CSV is saved in the working folder.
export_csv(all_orfs, DEFAULT_OUTPUT_FILE)
print(f"\nCSV saved to: {Path(DEFAULT_OUTPUT_FILE).resolve()}")



