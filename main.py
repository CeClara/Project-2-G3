
# ============================================================
# PROJECT 2: ORF FINDER TOOL
# Single-file version: core functions + Streamlit website
# Run with: streamlit run main.py
# ============================================================

import csv
import io
import statistics
from collections import Counter

import streamlit as st


# ============================================================
# SECTION 1: IMPORTS AND SHARED SETTINGS
# ============================================================

START_CODONS = ["ATG"]
STOP_CODONS = ["TAA", "TAG", "TGA"]
VALID_DNA_BASES = ["A", "C", "G", "T"]

MIN_ORF_LENGTH = 90
DEFAULT_OUTPUT_FILE = "orf_results.csv"

REVERSE_COMP_DNA = {
    "A": "T", "T": "A", "C": "G", "G": "C"
}

# Standard genetic code: DNA codon -> amino acid
GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# ============================================================
# SECTION 2: DNA INPUT AND VALIDATION
# ============================================================

def validate_dna(dnaseq):
    """Return cleaned uppercase DNA if valid, otherwise False."""

    if not isinstance(dnaseq, str):
        return False

    clean_seq = "".join(dnaseq.split()).upper()

    if not clean_seq:
        return False

    if any(base not in VALID_DNA_BASES for base in clean_seq):
        return False

    return clean_seq


def clean_file_text(file_text):
    """Remove FASTA headers if present, then validate the DNA."""

    sequence_lines = []

    for line in file_text.splitlines():
        line = line.strip()

        if not line or line.startswith(">"):
            continue

        sequence_lines.append(line)

    return validate_dna("".join(sequence_lines))


def load_dna(file_path):
    """Read a plain-text or FASTA DNA file and validate it."""

    try:
        with open(file_path, "r", encoding="utf-8") as file:
            file_text = file.read()

    except (OSError, UnicodeError) as error:
        st.error(f"Could not read the DNA file: {error}")
        return False

    return clean_file_text(file_text)


# ============================================================
# SECTION 3: REVERSE COMPLEMENT
# ============================================================

def reverse_complement(seq):
    """Return the reverse complement of a validated DNA sequence."""

    return "".join(
        REVERSE_COMP_DNA[base] for base in reversed(seq)
    )


# ============================================================
# SECTION 4: FIND ORFs IN ONE READING FRAME
# ============================================================

def find_orfs_in_frame(sequence, strand, frame, genome_length):
    """
    Find candidate ORFs in one reading frame.

    Coordinates are 1-based and inclusive.
    Reverse-strand coordinates refer to the original DNA:
    start is the lower coordinate and end is the higher.
    Length includes the stop codon.
    Nested starts before the same stop are reported separately.
    """

    orfs = []
    active_starts = []
    offset = frame - 1

    for pos in range(offset, len(sequence) - 2, 3):

        codon = sequence[pos:pos + 3]

        # Record a start codon.
        if codon in START_CODONS:
            active_starts.append(pos)

        # A stop codon closes active candidate ORFs.
        if codon in STOP_CODONS and active_starts:

            end_pos = pos + 2

            for start_pos in active_starts:

                dna_sequence = sequence[start_pos:end_pos + 1]

                # Convert to coordinates on the original DNA.
                if strand == "+":
                    start_coord = start_pos + 1
                    end_coord = end_pos + 1
                else:
                    start_coord = genome_length - end_pos
                    end_coord = genome_length - start_pos

                orfs.append({
                    "orf_id": "",
                    "start": start_coord,
                    "end": end_coord,
                    "strand": strand,
                    "frame": frame,
                    "length": len(dna_sequence),
                    "dna_sequence": dna_sequence,
                })

            active_starts = []

    return orfs


# ============================================================
# SECTION 5: SEARCH ALL SIX READING FRAMES
# ============================================================

def find_orfs(sequence):
    """Find candidate ORFs in all six reading frames."""

    sequence = validate_dna(sequence)

    if sequence is False:
        return []

    reverse_seq = reverse_complement(sequence)
    all_orfs = []

    for frame in (1, 2, 3):

        all_orfs.extend(
            find_orfs_in_frame(
                sequence, "+", frame, len(sequence)
            )
        )

        all_orfs.extend(
            find_orfs_in_frame(
                reverse_seq, "-", frame, len(sequence)
            )
        )

    # Sort by original DNA coordinates.
    all_orfs.sort(
        key=lambda orf: (
            orf["start"],
            orf["end"],
            orf["strand"],
            orf["frame"],
        )
    )

    # Assign unique IDs.
    for number, orf in enumerate(all_orfs, start=1):
        orf["orf_id"] = f"ORF_{number}"

    return all_orfs


# ============================================================
# SECTION 6: FILTER ORFs
# ============================================================

def filter_orfs(
    orfs,
    min_length=0,
    max_length=None,
    strand="both",
    frame="all",
):
    """Filter detected ORFs by length, strand and reading frame."""

    if strand not in ("both", "+", "-"):
        raise ValueError("strand must be 'both', '+' or '-'")

    if frame not in ("all", 1, 2, 3):
        raise ValueError("frame must be 'all', 1, 2 or 3")

    if min_length < 0:
        raise ValueError("min_length cannot be negative")

    if max_length is not None and max_length < min_length:
        raise ValueError(
            "max_length cannot be less than min_length"
        )

    filtered = []

    for orf in orfs:

        if orf["length"] < min_length:
            continue

        if max_length is not None and orf["length"] > max_length:
            continue

        if strand != "both" and orf["strand"] != strand:
            continue

        if frame != "all" and orf["frame"] != frame:
            continue

        filtered.append(orf)

    return filtered


# ============================================================
# SECTION 7: CALCULATE STATISTICS
# ============================================================

def calculate_statistics(orfs):
    """Calculate lengths, codon counts and amino acid counts."""

    lengths = [orf["length"] for orf in orfs]

    codon_counts = Counter()
    amino_acid_counts = Counter()

    for orf in orfs:

        dna = orf["dna_sequence"]

        # Exclude the terminal stop codon.
        coding_sequence = dna[:-3]

        for pos in range(0, len(coding_sequence) - 2, 3):

            codon = coding_sequence[pos:pos + 3]
            amino_acid = GENETIC_CODE.get(codon)

            if amino_acid is not None and amino_acid != "*":
                codon_counts[codon] += 1
                amino_acid_counts[amino_acid] += 1

    return {
        "number_of_orfs": len(orfs),
        "mean_length": (
            statistics.mean(lengths) if lengths else 0
        ),
        "shortest_length": min(lengths) if lengths else 0,
        "longest_length": max(lengths) if lengths else 0,
        "codon_counts": dict(codon_counts),
        "amino_acid_counts": dict(amino_acid_counts),
    }


# ============================================================
# SECTION 8: EXPORT RESULTS
# ============================================================

def export_results(orfs, file_path=DEFAULT_OUTPUT_FILE):
    """Write ORF information to a CSV file."""

    fieldnames = [
        "orf_id",
        "start",
        "end",
        "strand",
        "frame",
        "length",
        "dna_sequence",
    ]

    with open(
        file_path, "w", newline="", encoding="utf-8"
    ) as file:

        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(orfs)


def results_to_csv(orfs):
    """Create CSV data for the website's download button."""

    output = io.StringIO()

    fieldnames = [
        "orf_id",
        "start",
        "end",
        "strand",
        "frame",
        "length",
        "dna_sequence",
    ]

    writer = csv.DictWriter(output, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(orfs)

    return output.getvalue().encode("utf-8")


# ============================================================
# SECTION 9: STREAMLIT WEBSITE INTERFACE
# ============================================================

def main():

    st.set_page_config(
        page_title="ORF Finder",
        page_icon="🧬",
        layout="wide",
    )

    st.title("🧬 ORF Finder")
    st.write(
        "Find candidate open reading frames (ORFs) "
        "in a DNA sequence."
    )

    st.info(
        "This version uses ATG as the only start codon "
        "and TAA, TAG and TGA as stop codons."
    )

    input_method = st.radio(
        "Choose how to provide your DNA sequence:",
        ["Paste or type DNA", "Upload a DNA file"],
        horizontal=True,
    )

    # OPTION 1: Paste or type DNA.
    if input_method == "Paste or type DNA":

        dna_text = st.text_area(
            "DNA sequence",
            placeholder="Paste or type DNA using A, T, G and C...",
            height=160,
        )

        if st.button("Find ORFs", type="primary"):

            st.session_state.pop("orf_sequence", None)
            sequence = validate_dna(dna_text)

            if sequence is False:
                st.error(
                    "Enter a non-empty DNA sequence "
                    "containing only A, T, G and C."
                )
            else:
                st.session_state["orf_sequence"] = sequence

    # OPTION 2: Upload a DNA file.
    else:

        uploaded_file = st.file_uploader(
            "Drag and drop a DNA text or FASTA file, "
            "or browse to select one.",
            type=["txt", "fasta", "fa", "fna"],
        )

        if st.button("Find ORFs", type="primary"):

            st.session_state.pop("orf_sequence", None)

            if uploaded_file is None:
                st.error("Please upload a DNA file first.")

            else:
                try:
                    file_text = uploaded_file.getvalue().decode("utf-8")
                    sequence = clean_file_text(file_text)

                except UnicodeError:
                    sequence = False

                if sequence is False:
                    st.error(
                        "The file does not contain a valid DNA sequence. "
                        "Use A, T, G and C. FASTA headers are supported."
                    )
                else:
                    st.session_state["orf_sequence"] = sequence

    # Display results only after valid input has been submitted.
    if "orf_sequence" not in st.session_state:
        st.caption(
            "Enter a sequence or upload a file, then select Find ORFs."
        )
        return

    dna_sequence = st.session_state["orf_sequence"]

    st.success(
        f"DNA sequence loaded and validated "
        f"({len(dna_sequence):,} bases)."
    )

    all_orfs = find_orfs(dna_sequence)

    st.write(
        f"**Candidate ORFs detected before filtering:** "
        f"{len(all_orfs)}"
    )

    # Filtering controls.
    st.subheader("Filter ORFs")

    col1, col2, col3 = st.columns(3)

    with col1:
        min_length = st.number_input(
            "Minimum length (bases)",
            min_value=0,
            value=MIN_ORF_LENGTH,
            step=3,
        )

    with col2:
        max_length = st.number_input(
            "Maximum length (0 = no maximum)",
            min_value=0,
            value=0,
            step=3,
        )

    with col3:
        strand_label = st.selectbox(
            "Strand",
            ["Both", "Forward (+)", "Reverse (-)"],
        )

    frame_label = st.selectbox(
        "Reading frame",
        ["All", "Frame 1", "Frame 2", "Frame 3"],
    )

    strand_value = {
        "Both": "both",
        "Forward (+)": "+",
        "Reverse (-)": "-",
    }[strand_label]

    frame_value = (
        "all" if frame_label == "All"
        else int(frame_label[-1])
    )

    max_value = None if max_length == 0 else int(max_length)

    if max_value is not None and max_value < min_length:
        st.warning(
            "Maximum length must be greater than "
            "or equal to minimum length."
        )
        filtered_orfs = []

    else:
        filtered_orfs = filter_orfs(
            all_orfs,
            min_length=int(min_length),
            max_length=max_value,
            strand=strand_value,
            frame=frame_value,
        )

    # Summary statistics.
    stats = calculate_statistics(filtered_orfs)

    st.subheader("Summary statistics")

    metric1, metric2, metric3, metric4 = st.columns(4)

    metric1.metric(
        "ORFs after filtering",
        stats["number_of_orfs"],
    )

    metric2.metric(
        "Mean length",
        f'{stats["mean_length"]:.1f} bp',
    )

    metric3.metric(
        "Shortest length",
        f'{stats["shortest_length"]} bp',
    )

    metric4.metric(
        "Longest length",
        f'{stats["longest_length"]} bp',
    )

    # Results table and CSV download.
    st.subheader("ORF results")

    if filtered_orfs:

        st.dataframe(
            filtered_orfs,
            use_container_width=True,
            hide_index=True,
        )

        st.download_button(
            "Download results as CSV",
            data=results_to_csv(filtered_orfs),
            file_name=DEFAULT_OUTPUT_FILE,
            mime="text/csv",
        )

    else:
        st.info(
            "No ORFs match these filters. Try lowering "
            "the minimum length or changing the filters."
        )

    with st.expander("Codon counts"):
        st.json(stats["codon_counts"])

    with st.expander("Amino acid counts"):
        st.json(stats["amino_acid_counts"])


if __name__ == "__main__":
    main()
