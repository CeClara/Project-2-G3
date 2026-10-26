"""Streamlit web app for six-frame bacterial ORF finding."""
from __future__ import annotations

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from orf_finder.core import (clean_sequence, find_orfs, summary, codon_usage,
                             amino_acid_composition, intergenic_lengths)
from orf_finder.export import orfs_to_csv, orfs_to_fasta

st.set_page_config(page_title="Bacterial ORF Finder", layout="wide")

DISPLAY_COLS = ["id", "strand", "frame", "start", "end", "length_nt", "length_aa",
                "start_codon", "stop_codon", "nested", "overlap_count",
                "rbs_motif", "rbs_distance", "score", "confidence"]
COLORS = {"high": "#16803c", "medium": "#dd8a00", "low": "#c53b3b"}
MAX_MAP_ORFS = 5000

st.title("Bacterial ORF Finder")
st.write("Detect complete ORFs across all six reading frames, rank likely coding "
         "regions, and export results.")


# ---------------------------------------------------------------- analysis
def build_figures(orfs, seq_len):
    """Return the five Plotly figures for a list of ORFs."""
    empty = go.Figure()

    # Genome map: one trace per confidence level (much faster than 1 trace/ORF)
    map_fig = go.Figure()
    for conf, colour in COLORS.items():
        xs, ys, hover = [], [], []
        for o in orfs[:MAX_MAP_ORFS]:
            if o.confidence != conf:
                continue
            text = f"{o.id}: {o.start}-{o.end}, score {o.score}"
            xs += [o.start, o.end, None]
            ys += [o.frame, o.frame, None]
            hover += [text, text, None]
        if xs:
            map_fig.add_trace(go.Scatter(
                x=xs, y=ys, mode="lines", name=conf,
                line={"width": 8, "color": colour},
                hovertext=hover, hoverinfo="text"))
    map_fig.update_layout(title=f"ORF positions (first {MAX_MAP_ORFS:,} shown)",
                          xaxis_title="Genome coordinate", yaxis_title="Frame",
                          height=500)

    if orfs:
        lengths = pd.DataFrame({"length_nt": [o.length_nt for o in orfs],
                                "confidence": [o.confidence for o in orfs]})
        length_fig = px.histogram(lengths, x="length_nt", color="confidence",
                                  nbins=50, title="ORF length distribution")

        codons = pd.DataFrame(codon_usage(orfs))
        codon_fig = px.bar(codons.sort_values("frequency_per_1000", ascending=False),
                           x="codon", y="frequency_per_1000", color="amino_acid",
                           title="Codon usage")

        aa = pd.DataFrame(amino_acid_composition(orfs))
        aa_fig = px.bar(aa, x="amino_acid", y="percentage",
                        title="Amino-acid composition")

        gaps = intergenic_lengths(orfs, seq_len)
        gap_fig = px.histogram(x=gaps, nbins=40,
                               labels={"x": "Intergenic length (nt)"},
                               title="Intergenic length distribution")
    else:
        length_fig = codon_fig = aa_fig = gap_fig = empty

    return map_fig, length_fig, codon_fig, aa_fig, gap_fig


def run_analysis(text, min_nt, max_nt, strand, frames, starts, nested):
    seq = clean_sequence(text or "")
    orfs = find_orfs(seq, int(min_nt or 90), int(max_nt) if max_nt else None,
                     strand, set(frames) if frames else None,
                     set(starts or []), nested)
    records = [o.to_dict() for o in orfs]
    table = pd.DataFrame([{k: r[k] for k in DISPLAY_COLS} for r in records],
                         columns=DISPLAY_COLS)
    return {
        "seq_len": len(seq),
        "orfs": orfs,
        "table": table,
        "summary": summary(orfs, len(seq)),
        "figs": build_figures(orfs, len(seq)),
    }


# ------------------------------------------------------------------ inputs
left, right = st.columns(2)

with left:
    with st.container(border=True):
        st.subheader("1. Load DNA")
        uploaded = st.file_uploader("Drag a .txt/.fasta file here",
                                    type=["txt", "fasta", "fa", "fna"])
        # Copy a newly uploaded file into the text box (same behaviour as Dash)
        if uploaded is not None:
            upload_key = (uploaded.name, uploaded.size)
            if st.session_state.get("last_upload") != upload_key:
                st.session_state["last_upload"] = upload_key
                try:
                    st.session_state["sequence_input"] = uploaded.getvalue().decode("utf-8")
                except UnicodeDecodeError:
                    st.session_state["sequence_input"] = ""
                    st.error("Could not decode uploaded file. Please use a UTF-8 text file.")
        st.text_area("Or paste DNA/FASTA", key="sequence_input", height=130,
                     placeholder="ATG...TAA")

with right:
    with st.container(border=True):
        st.subheader("2. Detection settings")
        c1, c2 = st.columns(2)
        min_nt = c1.number_input("Minimum ORF length (nt)", value=90, min_value=6, step=3)
        max_nt = c2.number_input("Maximum length (blank = unlimited)", value=None,
                                 min_value=6, step=3)
        strand = st.radio("Strand", ["both", "+", "-"], horizontal=True,
                          format_func={"both": "Both", "+": "Forward", "-": "Reverse"}.get)
        frames = st.multiselect("Frames", [1, 2, 3, -1, -2, -3], placeholder="All frames")
        starts = st.multiselect("Start codons", ["ATG", "GTG", "TTG"],
                                default=["ATG", "GTG", "TTG"])
        nested = st.radio("Nested candidates", ["all", "longest"], horizontal=True,
                          format_func={"all": "Return all", "longest": "Longest per stop"}.get)
        run = st.button("Find ORFs", type="primary")

if run:
    try:
        with st.spinner("Searching for ORFs..."):
            st.session_state["result"] = run_analysis(
                st.session_state.get("sequence_input", ""),
                min_nt, max_nt, strand, frames, starts, nested)
        st.session_state.pop("error", None)
    except Exception as exc:  # noqa: BLE001
        st.session_state["result"] = None
        st.session_state["error"] = str(exc)

# ----------------------------------------------------------------- results
if st.session_state.get("error"):
    st.error(f"Error: {st.session_state['error']}")

res = st.session_state.get("result")
if res:
    orfs = res["orfs"]
    st.markdown(f"**Found {len(orfs):,} complete ORFs in {res['seq_len']:,} nt.**")

    # Summary cards
    items = list(res["summary"].items())
    per_row = 6
    for i in range(0, len(items), per_row):
        cols = st.columns(per_row)
        for col, (k, v) in zip(cols, items[i:i + per_row]):
            col.metric(k.replace("_", " ").title(), str(v))

    tab_table, tab_map, tab_len, tab_codon, tab_aa, tab_gap = st.tabs(
        ["ORF table", "Genome map", "Length distribution",
         "Codon usage", "Amino acids", "Intergenic lengths"])
    map_fig, length_fig, codon_fig, aa_fig, gap_fig = res["figs"]

    with tab_table:
        b1, b2, b3, _ = st.columns([1, 1, 1, 5])
        b1.download_button("Download CSV", orfs_to_csv(orfs),
                           "orf_results.csv", "text/csv")
        b2.download_button("DNA FASTA", orfs_to_fasta(orfs),
                           "orf_dna.fasta", "text/plain")
        b3.download_button("Protein FASTA", orfs_to_fasta(orfs, True),
                           "orf_proteins.fasta", "text/plain")

        table = res["table"]
        conf_filter = st.multiselect("Filter by confidence",
                                     ["high", "medium", "low"],
                                     default=["high", "medium", "low"])
        if not table.empty:
            table = table[table["confidence"].isin(conf_filter)]
        st.dataframe(table, use_container_width=True, hide_index=True)

    with tab_map:
        st.plotly_chart(map_fig, use_container_width=True)
    with tab_len:
        st.plotly_chart(length_fig, use_container_width=True)
    with tab_codon:
        st.plotly_chart(codon_fig, use_container_width=True)
    with tab_aa:
        st.plotly_chart(aa_fig, use_container_width=True)
    with tab_gap:
        st.plotly_chart(gap_fig, use_container_width=True)

st.caption("Scores are transparent heuristics for prioritisation, not proof of gene "
           "function. Validate predictions against an annotation or specialist gene predictor.")