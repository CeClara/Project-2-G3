# ============================================================
# PROJECT 2: ORF FINDER TOOL
# Run with: python -m streamlit run main2.py
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

VALID_DNA_BASES = set("ATGC")
MIN_ORF_LENGTH = 75
DEFAULT_OUTPUT_FILE = "orf_results.csv"

REVERSE_COMP_DNA = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
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
    """Remove whitespace, uppercase the sequence and validate it."""

    if not isinstance(dnaseq, str):
        return False

    sequence = "".join(dnaseq.split()).upper()

    if not sequence:
        return False

    if any(base not in VALID_DNA_BASES for base in sequence):
        return False

    return sequence


def parse_fasta(file_text):
    """Read plain DNA text or FASTA records separately."""

    lines = file_text.splitlines()
    has_headers = any(line.strip().startswith(">") for line in lines)

    if not has_headers:
        sequence = validate_dna(file_text)

        if sequence is False:
            return []

        return [("Sequence 1", sequence)]

    records = []
    current_name = None
    current_sequence = []

    for line in lines:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if current_name is not None:
                sequence = validate_dna("".join(current_sequence))

                if sequence is False:
                    raise ValueError(
                        f"Invalid or empty DNA sequence: {current_name}"
                    )

                records.append((current_name, sequence))

            current_name = line[1:].strip()

            if not current_name:
                current_name = f"Sequence {len(records) + 1}"

            current_sequence = []

        else:
            if current_name is None:
                raise ValueError(
                    "FASTA sequence is missing its header."
                )

            current_sequence.append(line)

    if current_name is not None:
        sequence = validate_dna("".join(current_sequence))

        if sequence is False:
            raise ValueError(
                f"Invalid or empty DNA sequence: {current_name}"
            )

        records.append((current_name, sequence))

    return records


def load_dna(file_path):
    """Load DNA records from a text or FASTA file."""

    try:
        with open(file_path, "r", encoding="utf-8") as file:
            return parse_fasta(file.read())

    except (OSError, UnicodeError, ValueError) as error:
        st.error(f"Could not read DNA file: {error}")
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


def translate_orf(dna_sequence, start_codons=None):
    """Translate DNA into amino acids, excluding a terminal stop."""

    if start_codons is None:
        start_codons = START_CODONS

    amino_acids = []

    for pos in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[pos:pos + 3]
        amino_acid = GENETIC_CODE.get(codon)

        if amino_acid is None or amino_acid == "*":
            break

        if pos == 0 and codon in start_codons:
            amino_acid = "M"

        amino_acids.append(amino_acid)

    return "".join(amino_acids)


# ============================================================
# SECTION 4: FIND ORFs IN ONE READING FRAME
# ============================================================

def find_orfs_in_frame(
    sequence,
    strand,
    frame,
    genome_length,
    start_codons=None,
    stop_codons=None,
    report_overlaps=True,
):
    """
    Find complete ORFs and partial ORFs with a start but no stop.

    report_overlaps=True:
        Report every start codon that reaches the same stop.

    report_overlaps=False:
        Keep only the longest ORF for each stop in this frame.

    Coordinates are 1-based and inclusive on the original sequence.
    Complete ORF lengths include the stop codon.
    Partial ORFs extend to the end of the available reading frame.
    """

    if start_codons is None:
        start_codons = START_CODONS

    if stop_codons is None:
        stop_codons = STOP_CODONS

    orfs = []
    offset = frame - 1
    active_starts = []

    aligned_end = offset + (
        (len(sequence) - offset) // 3
    ) * 3

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

        length = len(dna)

        if end_coord - start_coord + 1 != length:
            raise ValueError(
                f"Coordinate mismatch: {start_coord}-{end_coord}, "
                f"but sequence length is {length}."
            )

        orfs.append({
            "orf_id": "",
            "start": start_coord,
            "end": end_coord,
            "strand": strand,
            "frame": frame,
            "length": length,
            "gc_content": round(calculate_gc_content(dna), 2),
            "orf_type": orf_type,
            "dna_sequence": dna,
            "amino_acid_sequence": translate_orf(
                dna, start_codons
            ),
            "sequence_name": "",
        })

    for pos in range(offset, aligned_end, 3):
        codon = sequence[pos:pos + 3]

        if codon in start_codons:
            active_starts.append(pos)

        if codon in stop_codons and active_starts:
            if report_overlaps:
                starts_to_report = active_starts
            else:
                starts_to_report = [active_starts[0]]

            for start_pos in starts_to_report:
                add_orf(start_pos, pos + 3, "Complete")

            active_starts = []

    if active_starts:
        if report_overlaps:
            starts_to_report = active_starts
        else:
            starts_to_report = [active_starts[0]]

        for start_pos in starts_to_report:
            add_orf(start_pos, aligned_end, "Partial")

    return orfs


# ============================================================
# SECTION 5: SEARCH ALL SIX READING FRAMES
# ============================================================

def find_orfs(
    sequence,
    start_codons=None,
    stop_codons=None,
    report_overlaps=True,
):
    """Search all six reading frames."""

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
        all_orfs.extend(
            find_orfs_in_frame(
                sequence,
                "+",
                frame,
                len(sequence),
                start_codons,
                stop_codons,
                report_overlaps,
            )
        )

        all_orfs.extend(
            find_orfs_in_frame(
                reverse_sequence,
                "-",
                frame,
                len(sequence),
                start_codons,
                stop_codons,
                report_overlaps,
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
    """Filter ORFs by length, strand and reading frame."""

    if strand not in ("both", "+", "-"):
        raise ValueError("strand must be 'both', '+' or '-'")

    if frame not in ("all", 1, 2, 3):
        raise ValueError("frame must be 'all', 1, 2 or 3")

    if min_length < 0:
        raise ValueError("Minimum length cannot be negative")

    if max_length is not None and max_length < min_length:
        raise ValueError(
            "Maximum length cannot be less than minimum length"
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
    """Convert ORF results into downloadable CSV data."""

    output = io.StringIO()

    fieldnames = [
        "orf_id",
        "sequence_name",
        "start",
        "end",
        "strand",
        "frame",
        "length",
        "gc_content",
        "orf_type",
        "dna_sequence",
        "amino_acid_sequence",
    ]

    writer = csv.DictWriter(
        output,
        fieldnames=fieldnames,
        extrasaction="ignore",
    )

    writer.writeheader()
    writer.writerows(orfs)

    return output.getvalue().encode("utf-8")


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
        "Complete ORFs have a start codon followed by an in-frame "
        "stop codon. Partial ORFs have a start codon but no in-frame "
        "stop before the end of the available reading frame. "
        "Complete ORF lengths include the stop codon."
    )

    # DNA input
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
            sequence = validate_dna(dna_text)

            if sequence is False:
                st.error(
                    "Enter a non-empty DNA sequence containing "
                    "only A, T, G and C."
                )
                st.session_state.pop("orf_records", None)
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
            if uploaded_file is None:
                st.error("Please upload a DNA file first.")
                st.session_state.pop("orf_records", None)
            else:
                try:
                    file_text = uploaded_file.getvalue().decode("utf-8")
                    records = parse_fasta(file_text)

                    if not records:
                        st.error("No valid DNA sequences were found.")
                        st.session_state.pop("orf_records", None)
                    else:
                        st.session_state["orf_records"] = records

                except (UnicodeError, ValueError) as error:
                    st.error(f"Could not read DNA file: {error}")
                    st.session_state.pop("orf_records", None)

    if "orf_records" not in st.session_state:
        st.caption(
            "Enter a DNA sequence or upload a file, then select Find ORFs."
        )
        return

    records = st.session_state["orf_records"]
    total_bases = sum(len(sequence) for _, sequence in records)

    st.success(
        f"Loaded {len(records)} sequence(s), "
        f"containing {total_bases:,} bases in total."
    )

    # Sequence statistics
    st.subheader("Sequence statistics")

    overall_sequence = "".join(
        sequence for _, sequence in records
    )

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
            "FASTA records are analysed separately so ORFs cannot "
            "cross between different records."
        )

    # Codon settings
    st.subheader("Codon settings")

    codon_col1, codon_col2 = st.columns(2)

    with codon_col1:
        start_codons = st.multiselect(
            "Start codons",
            options=AVAILABLE_START_CODONS,
            default=START_CODONS,
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

    # Overlapping ORF reporting
    st.subheader("Overlapping ORF reporting")

    overlap_mode = st.radio(
        "How should overlapping ORFs be reported?",
        options=[
            "Keep all candidates (including overlaps)",
            "Keep only the longest ORF for each stop codon",
        ],
        index=0,
        help=(
            "Keep all candidates reports every start codon that reaches "
            "a stop. Longest-only keeps the earliest start and therefore "
            "the longest ORF for each stop codon within each reading frame. "
            "Different frames or strands can still overlap."
        ),
    )

    report_overlaps = (
        overlap_mode == "Keep all candidates (including overlaps)"
    )

    # Find ORFs in each record
    all_orfs = []

    for record_name, sequence in records:
        record_orfs = find_orfs(
            sequence,
            start_codons=start_codons,
            stop_codons=stop_codons,
            report_overlaps=report_overlaps,
        )

        for orf in record_orfs:
            orf["sequence_name"] = record_name

        all_orfs.extend(record_orfs)

    # Give ORFs unique IDs across all records
    for number, orf in enumerate(all_orfs, start=1):
        orf["orf_id"] = f"ORF_{number}"

    st.write(
        f"**Candidate ORFs detected before filtering:** {len(all_orfs)}"
    )

    complete_count = sum(
        orf["orf_type"] == "Complete" for orf in all_orfs
    )

    partial_count = sum(
        orf["orf_type"] == "Partial" for orf in all_orfs
    )

    count_col1, count_col2 = st.columns(2)
    count_col1.metric("Complete ORFs found", complete_count)
    count_col2.metric("Partial ORFs found", partial_count)

    # Filters
    st.subheader("Filter ORFs")

    orf_types = st.multiselect(
        "ORF type",
        options=["Complete", "Partial"],
        # FIX: Show partial ORFs by default as well.
        # This prevents boundary ORFs from being hidden.
        default=["Complete", "Partial"],
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
            "Maximum length must be greater than or equal "
            "to minimum length."
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

    # Summary statistics
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

    # Results table
    st.subheader("ORF results")

    if filtered_orfs:
        show_sequences = st.checkbox(
            "Show DNA and protein sequences in the table",
            value=True,
        )

        if show_sequences:
            table_orfs = filtered_orfs
        else:
            table_orfs = [
                {
                    key: value
                    for key, value in orf.items()
                    if key not in (
                        "dna_sequence",
                        "amino_acid_sequence",
                    )
                }
                for orf in filtered_orfs
            ]

        st.dataframe(
            table_orfs,
            use_container_width=True,
            hide_index=True,
        )

        st.download_button(
            "Download results as CSV",
            data=results_to_csv(filtered_orfs),
            file_name=DEFAULT_OUTPUT_FILE,
            mime="text/csv",
        )

        # Inspect an individual ORF
        st.subheader("Inspect an ORF")

        selected_id = st.selectbox(
            "Select an ORF",
            [orf["orf_id"] for orf in filtered_orfs],
        )

        selected_orf = next(
            orf for orf in filtered_orfs
            if orf["orf_id"] == selected_id
        )

        detail1, detail2, detail3, detail4 = st.columns(4)

        detail1.metric("Strand", selected_orf["strand"])
        detail2.metric("Frame", selected_orf["frame"])
        detail3.metric("Length", f'{selected_orf["length"]} bases')
        detail4.metric(
            "Coordinates",
            f'{selected_orf["start"]}–{selected_orf["end"]}',
        )

        st.write(f'**Type:** {selected_orf["orf_type"]}')

        st.write("**DNA sequence (coding strand orientation)**")
        st.code(selected_orf["dna_sequence"], language="text")

        st.write("**Protein translation**")
        st.code(
            selected_orf["amino_acid_sequence"] or "(empty)",
            language="text",
        )

    else:
        st.info(
            "No ORFs match these filters. Try selecting both ORF types, "
            "lowering the minimum length, or changing the strand or frame."
        )

    # Codon and amino-acid counts
    with st.expander("Codon counts"):
        st.json(stats["codon_counts"])

    with st.expander("Amino acid counts"):
        st.json(stats["amino_acid_counts"])


if __name__ == "__main__":
    main()