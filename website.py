

# ============================================================
# ORF FINDER - STREAMLIT WEBSITE INTERFACE
# Save this file as website.py in the same folder as main.py
# Run with: python -m streamlit run website.py
# ============================================================

import streamlit as st

# Import the ORF-finding functions and settings from main.py.
# Keep main.py in the same folder as website.py.
from main import (
    MIN_ORF_LENGTH,
    DEFAULT_OUTPUT_FILE,
    validate_dna,
    clean_file_text,
    find_orfs,
    filter_orfs,
    calculate_statistics,
    results_to_csv,
)


def main():
    st.set_page_config(
        page_title="ORF Finder",
        page_icon="🧬",
        layout="wide",
    )

    st.title("🧬 ORF Finder")
    st.write(
        "Find candidate open reading frames (ORFs) in a DNA sequence."
    )

    st.info(
        "This version uses ATG as the only start codon "
        "and TAA, TAG and TGA as stop codons."
    )

    # Choose how to provide the DNA sequence.
    input_method = st.radio(
        "Choose how to provide your DNA sequence:",
        ["Paste or type DNA", "Upload a DNA file"],
        horizontal=True,
    )

    # --------------------------------------------------------
    # OPTION 1: Paste or type DNA
    # --------------------------------------------------------
    if input_method == "Paste or type DNA":

        dna_text = st.text_area(
            "DNA sequence",
            placeholder="Paste or type a sequence using A, T, G and C...",
            height=160,
        )

        if st.button("Find ORFs", type="primary"):

            # Clear old results before checking new input.
            st.session_state.pop("orf_sequence", None)

            sequence = validate_dna(dna_text)

            if sequence is False:
                st.error(
                    "Please enter a non-empty DNA sequence "
                    "containing only A, T, G and C."
                )
            else:
                st.session_state["orf_sequence"] = sequence

    # --------------------------------------------------------
    # OPTION 2: Upload a DNA file
    # --------------------------------------------------------
    else:

        uploaded_file = st.file_uploader(
            "Drag and drop a DNA text or FASTA file here, "
            "or browse to select one",
            type=["txt", "fasta", "fa", "fna"],
        )

        if st.button("Find ORFs", type="primary"):

            # Clear old results before checking the new file.
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
                        "Use A, T, G and C; FASTA headers are supported."
                    )
                else:
                    st.session_state["orf_sequence"] = sequence

    # --------------------------------------------------------
    # Stop if no valid sequence has been submitted.
    # --------------------------------------------------------
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

    # --------------------------------------------------------
    # Find ORFs across all six reading frames.
    # --------------------------------------------------------
    all_orfs = find_orfs(dna_sequence)

    st.write(
        f"**Candidate ORFs detected before filtering:** {len(all_orfs)}"
    )

    # --------------------------------------------------------
    # Filters
    # --------------------------------------------------------
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

    # Convert website selections into values used by the functions.
    strand_value = {
        "Both": "both",
        "Forward (+)": "+",
        "Reverse (-)": "-",
    }[strand_label]

    if frame_label == "All":
        frame_value = "all"
    else:
        frame_value = int(frame_label[-1])

    if max_length == 0:
        max_value = None
    else:
        max_value = int(max_length)

    # Apply filters.
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

    # --------------------------------------------------------
    # Summary statistics
    # --------------------------------------------------------
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

    # --------------------------------------------------------
    # ORF results table and CSV download
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
            "No ORFs match these filters. "
            "Try lowering the minimum length or changing the filters."
        )

    # --------------------------------------------------------
    # Additional statistics
    # --------------------------------------------------------
    with st.expander("Codon counts"):
        st.json(stats["codon_counts"])

    with st.expander("Amino acid counts"):
        st.json(stats["amino_acid_counts"])


if __name__ == "__main__":
    main()
