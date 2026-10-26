"""Bacterial ORF Finder - single-file Streamlit app.

Run with:   python -m streamlit run orf_finder_app.py
Needs:      pip install streamlit pandas plotly

Sections
  1. ORF finding code   (your scanner, rewritten to be faster, plus scoring)
  2. Statistics         (codon usage, amino acids, intergenic gaps, summary)
  3. Export helpers     (CSV / FASTA)
  4. Plots
  5. Streamlit interface

Coordinates are 1-based, inclusive, on the forward strand of the input
(reverse-strand ORFs are converted back to forward coordinates).
"""
from __future__ import annotations

import csv
import io
import re
from bisect import bisect_left, bisect_right
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from typing import Optional

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# ======================================================================
# 1. ORF FINDING
# ======================================================================
DEFAULT_STARTS = ("ATG", "GTG", "TTG")
STOP_CODONS = ("TAA", "TAG", "TGA")

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {"".join(c): _AMINO[i] for i, c in enumerate(product(_BASES, repeat=3))}

# Shine-Dalgarno core motifs, strongest first
RBS_MOTIFS = ("AGGAGG", "GGAGG", "AGGAG", "GGAG", "GAGG", "AGGA")
_RBS_PATTERNS = tuple((m, re.compile(f"(?={m})")) for m in RBS_MOTIFS)

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

DISPLAY_COLS = ["id", "strand", "frame", "start", "end", "length_nt", "length_aa",
                "start_codon", "stop_codon", "nested", "overlap_count",
                "rbs_motif", "rbs_distance", "score", "confidence"]


@dataclass
class ORF:
    id: str
    strand: str           # "+" or "-"
    frame: int            # +1..+3 forward, -1..-3 reverse
    start: int            # 1-based, forward-strand coordinates
    end: int
    length_nt: int        # includes the stop codon
    length_aa: int        # excludes the stop codon
    start_codon: str
    stop_codon: str
    nested: bool          # True if an in-frame start inside a longer ORF
    overlap_count: int    # other ORFs overlapping this one (same-stop siblings excluded)
    rbs_motif: str
    rbs_distance: Optional[int]
    score: float
    confidence: str


def get_reverse_complement(sequence: str) -> str:
    """Reverse complement (str.translate is much faster than a per-base loop)."""
    return sequence.upper().translate(_COMPLEMENT)[::-1]


def clean_sequence(text: str) -> str:
    """Return an upper-case ACGTN string from raw DNA or single-record FASTA text."""
    lines = [ln.strip() for ln in (text or "").splitlines() if ln.strip()]
    if not lines:
        raise ValueError("No sequence provided.")
    seq_lines, records = [], 0
    for ln in lines:
        if ln.startswith(">"):
            records += 1
            if records > 1:          # only the first FASTA record is used
                break
            continue
        seq_lines.append(ln)
    seq = re.sub(r"\s+", "", "".join(seq_lines)).upper().replace("U", "T")
    bad = sorted(set(seq) - set("ACGTN"))
    if bad:
        raise ValueError(f"Invalid characters in sequence: {', '.join(bad[:5])}")
    if len(seq) < 6:
        raise ValueError("Sequence is too short.")
    return seq


def _translate(coding: str) -> str:
    """Translate a coding sequence (no stop codon). First residue is always M."""
    aa = [CODON_TABLE.get(coding[i:i + 3], "X") for i in range(0, len(coding) - 2, 3)]
    if aa:
        aa[0] = "M"      # GTG/TTG starts are read as Met in bacteria
    return "".join(aa)


def _codon_positions(seq: str, codons) -> tuple[list[int], list[int], list[int]]:
    """Positions of every codon in `codons`, bucketed by position % 3.

    Bucket n is exactly the set of codons in reading frame n, so one regex pass
    over the sequence replaces three slow Python loops over every codon.
    """
    pattern = re.compile("(?=(?:" + "|".join(sorted(codons)) + "))")
    buckets: tuple[list[int], list[int], list[int]] = ([], [], [])
    for m in pattern.finditer(seq):
        p = m.start()
        buckets[p % 3].append(p)
    return buckets


@lru_cache(maxsize=1 << 16)
def _rbs_in_window(upstream: str):
    """Shine-Dalgarno search in the (up to 20 nt) window before a start codon.

    Returns (motif, spacer) where spacer is the number of bases between the
    motif and the start codon (accepted range 4-14), or ("", None).
    """
    n = len(upstream)
    for motif, pattern in _RBS_PATTERNS:
        best = None
        for m in pattern.finditer(upstream):
            dist = n - (m.start() + len(motif))
            if 4 <= dist <= 14 and (best is None or abs(dist - 7) < abs(best - 7)):
                best = dist
        if best is not None:
            return motif, best
    return "", None


def _annotate_overlaps(orfs: list[ORF]) -> None:
    """overlap_count = other ORFs overlapping this one, ignoring ORFs that share its stop."""
    if not orfs:
        return
    starts = sorted(o.start for o in orfs)
    ends = sorted(o.end for o in orfs)
    stop_key = lambda o: (o.strand, o.end if o.strand == "+" else o.start)  # noqa: E731
    siblings = Counter(stop_key(o) for o in orfs)
    for o in orfs:
        overlapping = bisect_right(starts, o.end) - bisect_left(ends, o.start) - 1
        o.overlap_count = max(0, overlapping - (siblings[stop_key(o)] - 1))


def _score(o: ORF) -> float:
    """Transparent 0-100 heuristic: length, start codon, RBS, nesting, overlaps."""
    length = min(o.length_aa / 250, 1) * 45
    start = {"ATG": 15, "GTG": 7, "TTG": 4}.get(o.start_codon, 0)
    rbs = 0
    if o.rbs_motif:
        rbs = {6: 25, 5: 20, 4: 12}.get(len(o.rbs_motif), 0)
        if o.rbs_distance is not None and 5 <= o.rbs_distance <= 10:
            rbs += 5
    bonus = 0 if o.nested else 5
    penalty = min(o.overlap_count, 5) * 2
    return round(max(0.0, min(100.0, length + start + rbs + bonus - penalty)), 1)


def _confidence(score: float) -> str:
    return "high" if score >= 65 else "medium" if score >= 40 else "low"


def find_orfs(seq: str, min_nt: int = 90, max_nt: Optional[int] = None,
              strand: str = "both", frames=None, starts=None,
              nested: str = "all") -> list[ORF]:
    """Find complete ORFs (start codon ... stop codon) in all six frames.

    nested="all"      -> every in-frame start codon before a stop is reported
    nested="longest"  -> only the longest ORF for each stop codon
    """
    starts = set(DEFAULT_STARTS) if starts is None else set(starts)
    if not starts:
        raise ValueError("Select at least one start codon.")

    L = len(seq)
    orfs: list[ORF] = []
    for sym in ("+", "-"):
        if strand not in ("both", sym):
            continue
        s = seq if sym == "+" else get_reverse_complement(seq)
        stop_by_frame = _codon_positions(s, STOP_CODONS)
        start_by_frame = _codon_positions(s, starts)

        for offset in range(3):
            signed = offset + 1 if sym == "+" else -(offset + 1)
            if frames and signed not in frames:
                continue
            start_list = start_by_frame[offset]
            lo = 0
            # Walk the stop codons of this frame; every start codon found since
            # the previous stop belongs to the ORF that ends at this stop.
            for stop in stop_by_frame[offset]:
                hi = bisect_left(start_list, stop, lo)
                if hi > lo:
                    end0 = stop + 3
                    first = start_list[lo]
                    candidates = (first,) if nested == "longest" else start_list[lo:hi]
                    stop_codon = s[stop:end0]
                    for p in candidates:
                        length = end0 - p
                        if length < min_nt or (max_nt and length > max_nt):
                            continue
                        gstart, gend = (p + 1, end0) if sym == "+" else (L - end0 + 1, L - p)
                        motif, dist = _rbs_in_window(s[max(0, p - 20):p])
                        orfs.append(ORF(
                            id="", strand=sym, frame=signed, start=gstart, end=gend,
                            length_nt=length, length_aa=length // 3 - 1,
                            start_codon=s[p:p + 3], stop_codon=stop_codon,
                            nested=(p != first), overlap_count=0,
                            rbs_motif=motif, rbs_distance=dist,
                            score=0.0, confidence="low"))
                lo = hi

    orfs.sort(key=lambda o: (o.start, o.end, o.strand))
    _annotate_overlaps(orfs)
    for i, o in enumerate(orfs, 1):
        o.id = f"ORF_{i:05d}"
        o.score = _score(o)
        o.confidence = _confidence(o.score)
    return orfs


def orf_dna(o: ORF, seq: str) -> str:
    """Coding-strand DNA (5'->3') of an ORF, cut from the original sequence."""
    sub = seq[o.start - 1:o.end]
    return sub if o.strand == "+" else sub.translate(_COMPLEMENT)[::-1]


def orf_protein(o: ORF, seq: str) -> str:
    return _translate(orf_dna(o, seq)[:-3])


# ======================================================================
# 2. STATISTICS
# ======================================================================
def representatives(orfs: list[ORF]) -> list[ORF]:
    """One ORF per stop codon (the longest), so nested starts aren't double counted."""
    best: dict = {}
    for o in orfs:
        key = (o.strand, o.end if o.strand == "+" else o.start)
        if key not in best or o.length_nt > best[key].length_nt:
            best[key] = o
    return sorted(best.values(), key=lambda o: (o.start, o.end))


def _covered_bases(orfs: list[ORF]) -> int:
    total, cur_start, cur_end = 0, None, None
    for o in sorted(orfs, key=lambda o: o.start):
        if cur_end is None or o.start > cur_end + 1:
            if cur_end is not None:
                total += cur_end - cur_start + 1
            cur_start, cur_end = o.start, o.end
        else:
            cur_end = max(cur_end, o.end)
    if cur_end is not None:
        total += cur_end - cur_start + 1
    return total


def summary(orfs: list[ORF], seq_len: int) -> dict:
    n = len(orfs)
    conf = Counter(o.confidence for o in orfs)
    lengths = [o.length_nt for o in orfs]
    return {
        "total_orfs": n,
        "genome_length_nt": seq_len,
        "forward_orfs": sum(o.strand == "+" for o in orfs),
        "reverse_orfs": sum(o.strand == "-" for o in orfs),
        "high_confidence": conf["high"],
        "medium_confidence": conf["medium"],
        "low_confidence": conf["low"],
        "mean_length_nt": round(sum(lengths) / n, 1) if n else 0,
        "longest_orf_nt": max(lengths) if n else 0,
        "coding_density_pct": round(100 * _covered_bases(orfs) / seq_len, 1) if seq_len else 0,
    }


def composition(reps: list[ORF], seq: str) -> tuple[list[dict], list[dict]]:
    """Codon usage and amino-acid composition, from one pass over the ORFs."""
    codons: Counter = Counter()
    for o in reps:
        body = orf_dna(o, seq)[:-3]
        codons.update(map("".join, zip(body[0::3], body[1::3], body[2::3])))

    aa_counts: Counter = Counter()
    for codon, n in codons.items():
        aa_counts[CODON_TABLE.get(codon, "X")] += n
    for o in reps:                      # alternative starts are translated as Met
        if o.start_codon != "ATG":
            aa_counts[CODON_TABLE.get(o.start_codon, "X")] -= 1
            aa_counts["M"] += 1
    aa_counts.pop("X", None)

    total_c = sum(codons.values()) or 1
    total_a = sum(aa_counts.values()) or 1
    codon_rows = [{"codon": c, "amino_acid": CODON_TABLE.get(c, "X"), "count": n,
                   "frequency_per_1000": round(n / total_c * 1000, 2)}
                  for c, n in sorted(codons.items())]
    aa_rows = [{"amino_acid": a, "count": n, "percentage": round(n / total_a * 100, 2)}
               for a, n in sorted(aa_counts.items()) if n > 0]
    return codon_rows, aa_rows


def intergenic_lengths(reps: list[ORF]) -> list[int]:
    """Gap between consecutive ORFs (negative values mean the ORFs overlap)."""
    gaps, max_end = [], None
    for o in reps:
        if max_end is not None:
            gaps.append(o.start - max_end - 1)
        max_end = o.end if max_end is None else max(max_end, o.end)
    return gaps


# ======================================================================
# 3. EXPORT
# ======================================================================
def orfs_to_csv(orfs: list[ORF]) -> str:
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(DISPLAY_COLS)
    for o in orfs:
        d = vars(o)
        writer.writerow([d[c] for c in DISPLAY_COLS])
    return buf.getvalue()


def orfs_to_fasta(orfs: list[ORF], seq: str, protein: bool = False) -> str:
    lines = []
    for o in orfs:
        lines.append(f">{o.id} strand={o.strand} frame={o.frame} {o.start}-{o.end} "
                     f"len_nt={o.length_nt} score={o.score} confidence={o.confidence}")
        body = orf_protein(o, seq) if protein else orf_dna(o, seq)
        lines.extend(body[i:i + 60] for i in range(0, len(body), 60))
    return "\n".join(lines) + "\n" if lines else ""


# ======================================================================
# 4. PLOTS
# ======================================================================
COLORS = {"high": "#16803c", "medium": "#dd8a00", "low": "#c53b3b"}
MAX_MAP_ORFS = 5000


def build_figures(orfs, reps, codon_rows, aa_rows, gaps):
    empty = go.Figure()

    # Genome map: one trace per confidence level (much faster than one trace per ORF)
    map_fig = go.Figure()
    shown = reps[:MAX_MAP_ORFS]
    for conf, colour in COLORS.items():
        xs, ys, hover = [], [], []
        for o in shown:
            if o.confidence != conf:
                continue
            text = f"{o.id}: {o.start}-{o.end}, score {o.score}"
            xs += [o.start, o.end, None]
            ys += [o.frame, o.frame, None]
            hover += [text, text, None]
        if xs:
            map_fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name=conf,
                                         line={"width": 8, "color": colour},
                                         hovertext=hover, hoverinfo="text"))
    map_fig.update_layout(
        title=f"ORF positions (longest ORF per stop; first {MAX_MAP_ORFS:,} shown)",
        xaxis_title="Genome coordinate", yaxis_title="Frame", height=500)

    if not orfs:
        return map_fig, empty, empty, empty, empty

    lengths = pd.DataFrame({"length_nt": [o.length_nt for o in orfs],
                            "confidence": [o.confidence for o in orfs]})
    length_fig = px.histogram(lengths, x="length_nt", color="confidence", nbins=50,
                              color_discrete_map=COLORS, title="ORF length distribution")

    codons = pd.DataFrame(codon_rows).sort_values("frequency_per_1000", ascending=False)
    codon_fig = px.bar(codons, x="codon", y="frequency_per_1000", color="amino_acid",
                       title="Codon usage")

    aa_fig = px.bar(pd.DataFrame(aa_rows), x="amino_acid", y="percentage",
                    title="Amino-acid composition")

    gap_fig = (px.histogram(x=gaps, nbins=40, labels={"x": "Intergenic length (nt)"},
                            title="Intergenic length distribution")
               if gaps else empty)
    return map_fig, length_fig, codon_fig, aa_fig, gap_fig


def run_analysis(text, min_nt, max_nt, strand, frames, starts, nested) -> dict:
    seq = clean_sequence(text or "")
    orfs = find_orfs(seq, int(min_nt or 90), int(max_nt) if max_nt else None,
                     strand, set(frames) if frames else None,
                     set(starts) if starts else set(), nested)
    reps = representatives(orfs)
    codon_rows, aa_rows = composition(reps, seq) if orfs else ([], [])
    gaps = intergenic_lengths(reps)
    table = pd.DataFrame([vars(o) for o in orfs], columns=DISPLAY_COLS)
    return {
        "seq_len": len(seq),
        "orfs": orfs,
        "table": table,
        "summary": summary(orfs, len(seq)),
        "figs": build_figures(orfs, reps, codon_rows, aa_rows, gaps),
        # export files are built once here, not on every Streamlit rerun
        "csv": orfs_to_csv(orfs),
        "dna_fasta": orfs_to_fasta(orfs, seq),
        "protein_fasta": orfs_to_fasta(orfs, seq, protein=True),
    }


# ======================================================================
# 5. STREAMLIT INTERFACE
# ======================================================================
def main() -> None:
    st.set_page_config(page_title="Bacterial ORF Finder", layout="wide")
    st.title("Bacterial ORF Finder")
    st.write("Detect complete ORFs across all six reading frames, rank likely coding "
             "regions, and export results.")

    left, right = st.columns(2)

    with left:
        with st.container(border=True):
            st.subheader("1. Load DNA")
            uploaded = st.file_uploader("Drag a .txt/.fasta file here",
                                        type=["txt", "fasta", "fa", "fna"])
            # Copy a newly uploaded file into the text box
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
            starts = st.multiselect("Start codons", list(DEFAULT_STARTS),
                                    default=list(DEFAULT_STARTS))
            nested = st.radio("Nested candidates", ["all", "longest"], horizontal=True,
                              format_func={"all": "Return all",
                                           "longest": "Longest per stop"}.get)
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

    if st.session_state.get("error"):
        st.error(f"Error: {st.session_state['error']}")

    res = st.session_state.get("result")
    if res:
        orfs = res["orfs"]
        st.markdown(f"**Found {len(orfs):,} complete ORFs in {res['seq_len']:,} nt.**")

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
            b1.download_button("Download CSV", res["csv"], "orf_results.csv", "text/csv")
            b2.download_button("DNA FASTA", res["dna_fasta"], "orf_dna.fasta", "text/plain")
            b3.download_button("Protein FASTA", res["protein_fasta"],
                               "orf_proteins.fasta", "text/plain")

            table = res["table"]
            conf_filter = st.multiselect("Filter by confidence", ["high", "medium", "low"],
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


if __name__ == "__main__":
    main()