# the problem was the removal of ORF detection from our main code to this fixing it now 
# ============================================================
# PROJECT 2: ORF FINDER TOOL
# Run with: python -m streamlit run main.py
# ============================================================

import csv
import io
import statistics
from collections import Counter

import streamlit as st


# ============================================================
# SECTION 1: SETTINGS AND GENETIC CODE
# ============================================================

START_CODONS = ["ATG"]
STOP_CODONS = ["TAA", "TAG", "TGA"]

AVAILABLE_START_CODONS = ["ATG", "GTG", "TTG", "CTG"]
AVAILABLE_STOP_CODONS = ["TAA", "TAG", "TGA"]

VALID_DNA_BASES = ["A", "C", "G", "T"]

MIN_ORF_LENGTH = 90
DEFAULT_OUTPUT_FILE = "orf_results.csv"

REVERSE_COMP_DNA = {
    "A": "T", "T": "A", "C": "G", "G": "C"
}

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
    """Clean and validate a DNA sequence."""

    if not isinstance(dnaseq, str):
        return False

    sequence = "".join(dnaseq.split()).upper()

    if not sequence:
        return False

    if any(base not in VALID_DNA_BASES for base in sequence):
        return False

    return sequence


def parse_fasta(file_text):
    """
    Read plain DNA text or FASTA text.

    Each FASTA record is kept separate so ORFs cannot accidentally
    cross from one contig into another.
    """

    records = []
    current_name = "Sequence 1"
    current_sequence = []

    lines = file_text.splitlines()
    has_headers = any(line.strip().startswith(">") for line in lines)

    if not has_headers:
        sequence = validate_dna(file_text)
        if sequence is False:
            return []
        return [(current_name, sequence)]

    for line in lines:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if current_sequence:
                sequence = validate_dna("".join(current_sequence))
                if sequence is not False:
                    records.append((current_name, sequence))

            current_name = line[1:].strip() or f"Sequence {len(records) + 1}"
            current_sequence = []

        else:
            current_sequence.append(line)

    if current_sequence:
        sequence = validate_dna("".join(current_sequence))
        if sequence is not False:
            records.append((current_name, sequence))

    return records


def clean_file_text(file_text):
    """Return a single validated sequence from plain text or FASTA."""

    records = parse_fasta(file_text)

    if len(records) == 1:
        return records[0][1]

    return False


def load_dna(file_path):
    """Load and validate a DNA file."""

    try:
        with open(file_path, "r", encoding="utf-8") as file:
            return parse_fasta(file.read())

    except (OSError, UnicodeError) as error:
        st.error(f"Could not read the DNA file: {error}")
        return []


# ============================================================
# SECTION 3: SEQUENCE UTILITIES
# ============================================================

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""

    return "".join(
        REVERSE_COMP_DNA[base] for base in reversed(sequence)
    )


def calculate_gc_content(sequence):
    """Calculate the percentage of G and C bases."""

    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


# ============================================================
# SECTION 4: FIND ORFs IN ONE READING FRAME
# ============================================================
def translate_orf(dna_sequence, start_codons=None):
    """Translate complete codons, excluding the terminal stop.

    For alternative initiation codons, translate the first residue as M.
    Partial ORFs with no known start retain their ordinary first residue.
    """
    if start_codons is None:
        start_codons = START_CODONS
    amino_acids = []
    for pos in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[pos:pos + 3]
        amino_acid = GENETIC_CODE[codon]
        if amino_acid == "*":
            break
        if pos == 0 and codon in start_codons:
            amino_acid = "M"
        amino_acids.append(amino_acid)
    return "".join(amino_acids)


def find_orfs_in_frame(
    sequence,
    strand,
    frame,
    genome_length,
    start_codons=None,
    stop_codons=None,
):
    """Find complete, missing-start, missing-stop and missing-both ORFs.

    Frames are numbered 1..3 in the 5'-3' strand orientation.
    Coordinates are 1-based and inclusive in the original input DNA.
    Candidate partials are stop-delimited or extend to the sequence edge.
    """
    if start_codons is None:
        start_codons = START_CODONS
    if stop_codons is None:
        stop_codons = STOP_CODONS

    orfs = []
    offset = frame - 1
    # Only whole codons belong to the frame; exclude dangling 1-2 bases.
    aligned_end = offset + ((len(sequence) - offset) // 3) * 3
    if aligned_end <= offset:
        return orfs

    def add_orf(start_pos, end_exclusive, orf_type):
        dna = sequence[start_pos:end_exclusive]
        if not dna:
            return
        if strand == "+":
            start_coord = start_pos + 1
            end_coord = end_exclusive
        else:
            start_coord = genome_length - end_exclusive + 1
            end_coord = genome_length - start_pos
        orfs.append({
            "orf_id": "",
            "start": start_coord,
            "end": end_coord,
            "strand": strand,
            "frame": frame,
            "length": len(dna),
            "gc_content": round(calculate_gc_content(dna), 2),
            "orf_type": "Complete" if orf_type == "complete" else "Partial",
            "partial_type": orf_type,
            "dna_sequence": dna,
            "amino_acid_sequence": translate_orf(dna, start_codons),
            "nested": False,
        })

    active_starts = []
    segment_start = offset
    for pos in range(offset, aligned_end, 3):
        codon = sequence[pos:pos + 3]
        if codon in start_codons:
            active_starts.append(pos)
        if codon in stop_codons:
            if active_starts:
                for start_pos in active_starts:
                    add_orf(start_pos, pos + 3, "complete")
            elif pos > segment_start:
                # No start since last stop: candidate 5'-partial region.
                add_orf(segment_start, pos + 3, "missing_start")
            active_starts = []
            segment_start = pos + 3

    if active_starts:
        for start_pos in active_starts:
            add_orf(start_pos, aligned_end, "missing_stop")
    elif segment_start < aligned_end:
        add_orf(segment_start, aligned_end, "missing_both")

    # Nested = strictly contained in a longer candidate on same strand/frame.
    for orf in orfs:
        orf["nested"] = any(
            other is not orf
            and other["start"] <= orf["start"]
            and other["end"] >= orf["end"]
            and other["length"] > orf["length"]
            for other in orfs
        )
    return orfs


# ============================================================
# SECTION 5: SEARCH ALL SIX READING FRAMES
# ============================================================

def find_orfs(sequence, start_codons=None, stop_codons=None):
    """Find candidate ORFs on both strands in all three frames."""

    sequence = validate_dna(sequence)

    if sequence is False:
        return []

    if start_codons is None:
        start_codons = START_CODONS

    if stop_codons is None:
        stop_codons = STOP_CODONS

    reverse_sequence = reverse_complement(sequence)
    all_orfs = []

    for frame in (1, 2, 3):

        # Forward strand.
        all_orfs.extend(
            find_orfs_in_frame(
                sequence, "+", frame, len(sequence),
                start_codons, stop_codons
            )
        )

        # Reverse strand.
        all_orfs.extend(
            find_orfs_in_frame(
                reverse_sequence, "-", frame, len(sequence),
                start_codons, stop_codons
            )
        )

    all_orfs.sort(
        key=lambda orf: (
            orf["start"],
            orf["end"],
            orf["strand"],
            orf["frame"],
        )
    )

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
    """Filter ORFs by length, strand and frame."""

    if strand not in ("both", "+", "-"):
        raise ValueError("strand must be 'both', '+' or '-'")

    if frame not in ("all", 1, 2, 3):
        raise ValueError("frame must be 'all', 1, 2 or 3")

    if min_length < 0:
        raise ValueError("min_length cannot be negative")

    if max_length is not None and max_length < min_length:
        raise ValueError("max_length cannot be less than min_length")

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
# SECTION 7: STATISTICS
# ============================================================

def calculate_statistics(orfs):
    """Calculate ORF length, codon and amino-acid statistics."""

    lengths = [orf["length"] for orf in orfs]
    codon_counts = Counter()
    amino_acid_counts = Counter()

    for orf in orfs:

        dna = orf["dna_sequence"]

        if orf["orf_type"] == "Complete":
            coding_sequence = dna[:-3]
        else:
            coding_sequence = dna

        for pos in range(0, len(coding_sequence) - 2, 3):

            codon = coding_sequence[pos:pos + 3]
            amino_acid = GENETIC_CODE.get(codon)

            if amino_acid is not None and amino_acid != "*":
                codon_counts[codon] += 1
                amino_acid_counts[amino_acid] += 1

    return {
        "number_of_orfs": len(orfs),
        "mean_length": statistics.mean(lengths) if lengths else 0,
        "shortest_length": min(lengths) if lengths else 0,
        "longest_length": max(lengths) if lengths else 0,
        "codon_counts": dict(codon_counts),
        "amino_acid_counts": dict(amino_acid_counts),
    }


# ============================================================
# SECTION 8: CSV EXPORT
# ============================================================

def results_to_csv(orfs):
    """Convert ORF results to CSV data."""

    output = io.StringIO()

    fieldnames = [
        "orf_id",
        "start",
        "end",
        "strand",
        "frame",
        "length",
        "gc_content",
        "orf_type",
        "dna_sequence",
    ]

    writer = csv.DictWriter(output, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(orfs)

    return output.getvalue().encode("utf-8")


def export_results(orfs, file_path=DEFAULT_OUTPUT_FILE):
    """Save ORF results to a CSV file."""

    with open(file_path, "wb") as file:
        file.write(results_to_csv(orfs))


# ============================================================
# SECTION 9: STREAMLIT WEBSITE
# ============================================================

def main():

    st.set_page_config(
        page_title="ORF Finder",
        page_icon="🧬",
        layout="wide",
    )

    st.title("🧬 ORF Finder")

    st.write(
        "Find candidate open reading frames (ORFs) in DNA sequences."
    )

    st.info(
        "Complete ORFs have a start codon and an in-frame stop codon. "
        "Partial ORFs have a start codon without an in-frame stop "
        "before the end of that sequence or contig."
    )

    # --------------------------------------------------------
    # DNA input
    # --------------------------------------------------------

    input_method = st.radio(
        "Choose how to provide your DNA sequence:",
        ["Paste or type DNA", "Upload a DNA file"],
        horizontal=True,
    )

    if input_method == "Paste or type DNA":

        dna_text = st.text_area(
            "DNA sequence",
            placeholder="Paste or type DNA using A, T, G and C...",
            height=160,
        )

        if st.button("Find ORFs", type="primary"):

            st.session_state.pop("orf_records", None)

            sequence = validate_dna(dna_text)

            if sequence is False:
                st.error(
                    "Enter a non-empty DNA sequence containing "
                    "only A, T, G and C."
                )
            else:
                st.session_state["orf_records"] = [
                    ("Sequence 1", sequence)
                ]

    else:

        uploaded_file = st.file_uploader(
            "Upload a DNA text or FASTA file",
            type=["txt", "fasta", "fa", "fna"],
        )

        if st.button("Find ORFs", type="primary"):

            st.session_state.pop("orf_records", None)

            if uploaded_file is None:
                st.error("Please upload a DNA file first.")

            else:
                try:
                    file_text = uploaded_file.getvalue().decode("utf-8")
                    records = parse_fasta(file_text)

                except UnicodeError:
                    records = []

                if not records:
                    st.error(
                        "No valid DNA sequences were found. "
                        "Use A, T, G and C. FASTA headers are supported."
                    )
                else:
                    st.session_state["orf_records"] = records

    if "orf_records" not in st.session_state:
        st.caption(
            "Enter a DNA sequence or upload a file, then select Find ORFs."
        )
        return

    records = st.session_state["orf_records"]

    total_bases = sum(len(sequence) for _, sequence in records)

    st.success(
        f"Loaded {len(records)} sequence(s), containing "
        f"{total_bases:,} bases in total."
    )

    # --------------------------------------------------------
    # Sequence statistics
    # --------------------------------------------------------

    st.subheader("Sequence statistics")

    overall_sequence = "".join(sequence for _, sequence in records)

    gc_col, length_col = st.columns(2)

    gc_col.metric(
        "Overall GC content",
        f"{calculate_gc_content(overall_sequence):.2f}%",
    )

    length_col.metric(
        "Total sequence length",
        f"{total_bases:,} bases",
    )

    if len(records) > 1:
        st.caption(
            "FASTA records are analysed separately to prevent ORFs "
            "from crossing contig boundaries."
        )

    # --------------------------------------------------------
    # Custom start and stop codons
    # --------------------------------------------------------

    st.subheader("Codon settings")

    codon_col1, codon_col2 = st.columns(2)

    with codon_col1:
        start_codons = st.multiselect(
            "Start codons",
            options=AVAILABLE_START_CODONS,
            default=START_CODONS,
            help="ATG is the default. GTG, TTG and CTG are also available.",
        )

    with codon_col2:
        stop_codons = st.multiselect(
            "Stop codons",
            options=AVAILABLE_STOP_CODONS,
            default=STOP_CODONS,
        )

    if not start_codons:
        st.warning("Select at least one start codon.")
        return

    if not stop_codons:
        st.warning("Select at least one stop codon.")
        return

    # --------------------------------------------------------
    # Find ORFs in every sequence record separately
    # --------------------------------------------------------

    all_orfs = []

    for record_name, sequence in records:

        record_orfs = find_orfs(
            sequence,
            start_codons=start_codons,
            stop_codons=stop_codons,
        )

        for orf in record_orfs:
            orf["sequence_name"] = record_name

        all_orfs.extend(record_orfs)

    # Give every ORF a unique ID across all records.
    for number, orf in enumerate(all_orfs, start=1):
        orf["orf_id"] = f"ORF_{number}"

    st.write(
        f"**Candidate ORFs detected before filtering:** {len(all_orfs)}"
    )

    partial_count = sum(
        orf["orf_type"] == "Partial" for orf in all_orfs
    )

    complete_count = sum(
        orf["orf_type"] == "Complete" for orf in all_orfs
    )

    count_col1, count_col2 = st.columns(2)

    count_col1.metric("Complete ORFs found", complete_count)
    count_col2.metric("Partial ORFs found", partial_count)

    if not any(
        codon in sequence
        for _, sequence in records
        for codon in start_codons
    ) and not any(
        codon in reverse_complement(sequence)
        for _, sequence in records
        for codon in start_codons
    ):
        st.warning("No selected start codons were found in the sequences.")

    if not any(
        codon in sequence
        for _, sequence in records
        for codon in stop_codons
    ) and not any(
        codon in reverse_complement(sequence)
        for _, sequence in records
        for codon in stop_codons
    ):
        st.warning(
            "No selected stop codons were found. "
            "Complete ORFs may not be detected."
        )

    # --------------------------------------------------------
    # Filters
    # --------------------------------------------------------

    st.subheader("Filter ORFs")

    orf_types = st.multiselect(
        "ORF type",
        options=["Complete", "Partial"],
        default=["Complete"],
        help="Choose Complete, Partial, or both.",
    )

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
            "Maximum length must be greater than or equal to minimum length."
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

        filtered_orfs = [
            orf for orf in filtered_orfs
            if orf["orf_type"] in orf_types
        ]

    # --------------------------------------------------------
    # Summary statistics
    # --------------------------------------------------------

    stats = calculate_statistics(filtered_orfs)

    st.subheader("Summary statistics")

    metric1, metric2, metric3, metric4 = st.columns(4)

    metric1.metric("ORFs after filtering", stats["number_of_orfs"])

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

    # --------------------------------------------------------
    # Results table and CSV download
    # --------------------------------------------------------

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
            "No ORFs match these filters. Try selecting Complete, "
            "Partial, or both; lowering the minimum length; or "
            "changing the codon, strand or frame settings."
        )

    # --------------------------------------------------------
    # Codon and amino-acid counts
    # --------------------------------------------------------


    with st.expander("Codon counts"):
        st.json(stats["codon_counts"])

    with st.expander("Amino acid counts"):
        st.json(stats["amino_acid_counts"])

    if filtered_orfs:
        st.subheader("Selected ORF amino acid translation")
        chosen_id = st.selectbox("Select ORF", [o["orf_id"] for o in filtered_orfs])
        chosen = next(o for o in filtered_orfs if o["orf_id"] == chosen_id)
        st.code(chosen["amino_acid_sequence"], language="text")
        st.caption(f"Type: {chosen['partial_type']}; nested: {chosen['nested']}")


if __name__ == "__main__":
    main()
