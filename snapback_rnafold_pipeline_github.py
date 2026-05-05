#!/usr/bin/env python3
"""
Snapback breakpoint analysis with RNAfold. 

This script is adapted from my notebook so I can run the full workflow
from the command line and reuse it easily for multiple samples.

To run:
    python snapback_rnafold_pipeline_server.py

Make sure RNAfold (ViennaRNA) is installed and available in PATH.
You can check this using:
    which RNAfold

Also check all the path to files in File reference pair before running the code..
For additional samples, just add more entries in FILE_REFERENCE_PAIRS.
"""

import argparse
import re
import shutil
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import chisquare
import matplotlib.pyplot as plt


# ---- mention the path ----
PROJECT_DIR = Path(".")

FILE_REFERENCE_PAIRS = [
    (
        PROJECT_DIR / "example.snapback.tile.zmw.counts",
        PROJECT_DIR / "reference_split.fa",
        PROJECT_DIR / "sample_output_new",
    ),
]

# find RNAfold from PATH instead of hardcoding a server path
RNAFOLD_CMD = shutil.which("RNAfold")

if RNAFOLD_CMD is None:
    raise RuntimeError(
        "RNAfold not found in PATH.\n"
        "Please install ViennaRNA or add RNAfold to PATH.\n"
        "Check with: which RNAfold"
    )

RESET_OUTPUT_FOLDER = True  # True will recreate sample_output_new each time


WINDOWS: list[tuple[int, int, str]] = [
    (10, 10, "10bp-Brkp-10bp"),
    (15, 15, "15bp-Brkp-15bp"),
    (20, 20, "20bp-Brkp-20bp"),
    (0, 40, "40bp-Brkp"),
    (40, 0, "Brkp-40bp"),
    (5, 35, "5bp-Brkp-35bp"),
    (35, 5, "35bp-Brkp-5bp"),
    (10, 30, "10bp-Brkp-30bp"),
    (15, 25, "15bp-Brkp-25bp"),
    (25, 15, "25bp-Brkp-15bp"),
    (30, 10, "30bp-Brkp-10bp"),
]



# ---------------------------------------------------------
# STEP 1: Parse snapback tile-count file
# ---------------------------------------------------------
# This reads the subparser output file (.tile.zmw.counts).
# Each line has a pattern, count, and proportion. Here I pull
# out the two payload coordinates, strand information, counts,
# and proportions into a cleaner CSV file for later steps.
# ---------------------------------------------------------
def parse_snapback_tile_file(tile_file: Path, output_folder: Path) -> Path:
    # make a simple csv from the subparser tile count file
    output_csv = output_folder / f"{tile_file.stem}.csv"
    pattern = re.compile(r"Payload\S*\[(\d+)-(\d+)\]\(([tf\+\-])\)")

    rows = []
    with tile_file.open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip().replace("e-", "e")
            if not line:
                continue

            matches = pattern.findall(line)
            if len(matches) < 2:
                continue

            count, proportion = line.split()[0], line.split()[1]
            start_a, end_a, strand_a = matches[0]
            start_b, end_b, strand_b = matches[1]

            strand_a = "t" if strand_a in {"t", "+"} else "f"
            strand_b = "t" if strand_b in {"t", "+"} else "f"

            rows.append(
                {
                    "Count": float(count),
                    "Proportion": float(proportion),
                    "StartA": int(start_a),
                    "EndA": int(end_a),
                    "StrandA": strand_a,
                    "StartB": int(start_b),
                    "EndB": int(end_b),
                    "StrandB": strand_b,
                }
            )

    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    print(f"done: Parsed tile file: {output_csv}")
    return output_csv



# ---------------------------------------------------------
# STEP 2: Summarize breakpoint positions
# ---------------------------------------------------------
# After parsing the tile file, this step collects breakpoint
# positions separately for plus and minus strand. For this
# snapback output, I use EndA for the plus side and StartB
# for the minus side, then calculate count and proportion.
# ---------------------------------------------------------
def summarize_breakpoints(parsed_csv: Path, output_folder: Path) -> tuple[Path, Path]:
    # count breakpoints separately for plus and minus side
    df = pd.read_csv(parsed_csv)

    plus_df = (
        df.loc[df["StrandA"] == "t", ["Count", "EndA"]]
        .groupby("EndA", as_index=False)["Count"]
        .sum()
    )
    plus_total = plus_df["Count"].sum()
    plus_df["Proportion"] = plus_df["Count"] / plus_total if plus_total else 0
    plus_path = output_folder / "PROPORTION_enda.csv"
    plus_df.to_csv(plus_path, index=False)

    minus_df = (
        df.loc[df["StrandB"] == "t", ["Count", "StartB"]]
        .groupby("StartB", as_index=False)["Count"]
        .sum()
    )
    minus_total = minus_df["Count"].sum()
    minus_df["Proportion"] = minus_df["Count"] / minus_total if minus_total else 0
    minus_path = output_folder / "PROPORTION_startb.csv"
    minus_df.to_csv(minus_path, index=False)

    print(f"done: Wrote plus breakpoint table:  {plus_path}")
    print(f"done: Wrote minus breakpoint table: {minus_path}")
    return plus_path, minus_path



# ---------------------------------------------------------
# STEP 3A: Read payload sequence from reference FASTA
# ---------------------------------------------------------
# The reference FASTA should contain a sequence named Payload.
# If only one sequence is present, the script uses that sequence.
# This keeps the script flexible for small GitHub examples.
# ---------------------------------------------------------
def read_payload_from_reference(reference_fasta: Path) -> SeqRecord:
    # get payload sequence from reference fasta
    records = list(SeqIO.parse(reference_fasta, "fasta"))
    for record in records:
        if record.id == "Payload" or record.name == "Payload":
            return record

    if len(records) == 1:
        print("note: No record named 'Payload' found. Using the only FASTA record.")
        return records[0]

    names = ", ".join(record.id for record in records)
    raise ValueError(
        f"Could not find a FASTA record named 'Payload' in {reference_fasta}. "
        f"Available records: {names}"
    )



# ---------------------------------------------------------
# STEP 3B: Create plus and minus payload FASTA files
# ---------------------------------------------------------
# The plus strand is the original payload sequence.
# The minus strand is the reverse complement. These two FASTA
# files are used to extract local breakpoint sequence windows.
# ---------------------------------------------------------
def write_payload_fastas(reference_fasta: Path, output_folder: Path) -> tuple[Path, Path]:
    # write plus payload and reverse-complement sequence
    payload = read_payload_from_reference(reference_fasta)
    plus_record = SeqRecord(payload.seq, id="Payload", description="")
    minus_record = SeqRecord(
        payload.seq.reverse_complement(),
        id="Payload",
        description="Reverse complement",
    )

    plus_path = output_folder / "payload_plusstrand.fasta"
    minus_path = output_folder / "payload_minusstrand.fasta"
    SeqIO.write(plus_record, plus_path, "fasta")
    SeqIO.write(minus_record, minus_path, "fasta")

    print(f"done: Wrote plus payload FASTA:  {plus_path}")
    print(f"done: Wrote minus payload FASTA: {minus_path}")
    return plus_path, minus_path



# ---------------------------------------------------------
# Helper: read FASTA sequence
# ---------------------------------------------------------
# This small helper reads a FASTA file and returns the sequence
# as uppercase text so it can be indexed by coordinate.
# ---------------------------------------------------------
def read_fasta_sequence(fasta_path: Path) -> str:
    # read fasta as one sequence
    return "".join(str(record.seq).upper() for record in SeqIO.parse(fasta_path, "fasta"))



# ---------------------------------------------------------
# Helper: get sequence around one breakpoint
# ---------------------------------------------------------
# Coordinates from the tile file are 1-based. Python indexing is
# 0-based, so this function converts the coordinate and extracts
# the requested left/right sequence window around the breakpoint.
# ---------------------------------------------------------
def get_window_sequence(seq: str, breakpoint_1based: int, left: int, right: int) -> str:
    # breakpoint coordinate is 1-based. include the breakpoint base.
    idx = int(breakpoint_1based) - 1
    start = max(0, idx - left)
    end = min(len(seq), idx + right + 1)
    return seq[start:end]



# ---------------------------------------------------------
# STEP 4: Generate breakpoint-centered sequence windows
# ---------------------------------------------------------
# For each breakpoint, this writes multiple windows around that
# coordinate. I use symmetric and asymmetric windows because the
# strongest local RNAfold signal may not always be centered exactly
# the same way around the breakpoint.
# ---------------------------------------------------------
def write_breakpoint_windows(
    sequence_fasta: Path,
    breakpoint_csv: Path,
    breakpoint_column: str,
    output_fasta: Path,
) -> Path:
    # write fasta sequences around each breakpoint
    seq = read_fasta_sequence(sequence_fasta)
    df = pd.read_csv(breakpoint_csv)
    breakpoints = df[breakpoint_column].astype(int).tolist()

    with output_fasta.open("w") as out:
        for bp in breakpoints:
            for left, right, label in WINDOWS:
                window_seq = get_window_sequence(seq, bp, left, right)
                if not window_seq:
                    continue
                out.write(f"> {bp}_{label}\n{window_seq}\n")

    print(f"done: Wrote breakpoint windows: {output_fasta}")
    return output_fasta



# ---------------------------------------------------------
# STEP 5: Run RNAfold
# ---------------------------------------------------------
# RNAfold is called from Python using subprocess. It predicts
# local secondary structures and reports minimum free energy (MFE)
# for each breakpoint-window sequence.
# ---------------------------------------------------------
def run_rnafold(input_fasta: Path, output_file: Path, rnafold_cmd: str = RNAFOLD_CMD) -> Path:
    # RNAfold creates .ps files when -p is used.
    # I run RNAfold inside a small side folder so those files do not go into the main notebook folder.
    output_file.parent.mkdir(parents=True, exist_ok=True)
    ps_folder = output_file.parent / "rnafold_ps_files"
    ps_folder.mkdir(parents=True, exist_ok=True)

    with input_fasta.open("r") as infile:
        result = subprocess.run(
            [rnafold_cmd, "-p", "-d2", "--noLP"],
            stdin=infile,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            cwd=ps_folder,
        )

    output_file.write_text(result.stdout)

    if result.returncode != 0:
        raise RuntimeError(
            f"RNAfold failed for {input_fasta}.\n"
            f"Command: {rnafold_cmd} -p -d2 --noLP\n"
            f"Error:\n{result.stderr}"
        )

    print(f"done: RNAfold output: {output_file}")
    print(f"done: RNAfold .ps files: {ps_folder}")
    return output_file



# ---------------------------------------------------------
# STEP 6: Extract MFE values from RNAfold output
# ---------------------------------------------------------
# RNAfold output contains sequence, structure, and energy lines.
# This step extracts the MFE for each breakpoint window and keeps
# the strongest value per breakpoint for downstream comparison.
# ---------------------------------------------------------
def extract_mfe_from_rnafold(
    rnafold_file: Path,
    breakpoint_column: str,
    output_csv: Path,
) -> Path:
    # parse MFE values from RNAfold output. Keep strongest MFE per breakpoint.
    lines = rnafold_file.read_text().splitlines()
    rows = []

    for i in range(2, len(lines), 6):
        try:
            structure_1 = lines[i].split(" ")[0]
            structure_2 = lines[i + 2].split(" ")[0]
            header = lines[i - 2].strip()
            mfe_line = lines[i].strip()
        except IndexError:
            continue

        if structure_1 != structure_2:
            continue

        bp_match = re.search(r">\s*(\d+)_", header)
        mfe_match = re.search(r"\((-?\d+(?:\.\d+)?)\)", mfe_line)
        if not bp_match or not mfe_match:
            continue

        rows.append(
            {
                breakpoint_column: int(bp_match.group(1)),
                "MFE": abs(float(mfe_match.group(1))),
                "Sliding_window": header.replace(">", "").strip(),
            }
        )

    df = pd.DataFrame(rows)
    if df.empty:
        print(f"note: No MFE rows extracted from {rnafold_file}")
        df = pd.DataFrame(columns=[breakpoint_column, "MFE", "ranked", "Sliding_window"])
    else:
        df["ranked"] = df.groupby(breakpoint_column)["MFE"].rank("dense", ascending=False)
        df = df.loc[df["ranked"] == 1, [breakpoint_column, "MFE", "ranked", "Sliding_window"]]
        df = df.drop_duplicates(subset=breakpoint_column, keep="first")

    df.to_csv(output_csv, index=False)
    print(f"done: Wrote MFE table: {output_csv}")
    return output_csv



# ---------------------------------------------------------
# STEP 7: Merge breakpoint abundance with MFE
# ---------------------------------------------------------
# This combines breakpoint counts/proportions with the RNAfold MFE
# table for plus and minus strand. The final combined table is used
# for plotting and nucleotide analysis.
# ---------------------------------------------------------
def merge_plus_minus_mfe(output_folder: Path) -> Path:
    # combine count/proportion and mfe information
    plus_counts = pd.read_csv(output_folder / "PROPORTION_enda.csv")
    plus_mfe = pd.read_csv(output_folder / "Drp_duplicates_column_of_Intrst_plus.csv")
    minus_counts = pd.read_csv(output_folder / "PROPORTION_startb.csv")
    minus_mfe = pd.read_csv(output_folder / "Drp_duplicates_column_of_Intrst_minus.csv")

    plus_df = plus_counts.merge(plus_mfe[["EndA", "MFE"]], on="EndA", how="left")
    plus_df["MFE"] = plus_df["MFE"].fillna(0)
    plus_df["Category"] = "Coding / Plus strand"
    plus_df = plus_df.rename(columns={"EndA": "Breakpoint"})[
        ["Breakpoint", "Count", "Proportion", "MFE", "Category"]
    ]

    minus_df = minus_counts.merge(minus_mfe[["StartB", "MFE"]], on="StartB", how="left")
    minus_df["MFE"] = minus_df["MFE"].fillna(0)
    minus_df["Category"] = "Non-Coding / Minus strand"
    minus_df = minus_df.rename(columns={"StartB": "Breakpoint"})[
        ["Breakpoint", "Count", "Proportion", "MFE", "Category"]
    ]

    combined = pd.concat([plus_df, minus_df], ignore_index=True)
    out_path = output_folder / "plus_minus_MFE.csv"
    combined.to_csv(out_path, index=False)
    print(f"done: Wrote combined MFE table: {out_path}")
    return out_path



# ---------------------------------------------------------
# Helper: get nucleotide at breakpoint coordinate
# ---------------------------------------------------------
# This returns the base present at a 1-based breakpoint position.
# If the coordinate is outside the sequence, it returns N.
# ---------------------------------------------------------
def get_base_1based(seq: str, pos: int) -> str:
    # fasta coordinates are 1-based here
    p = int(pos)
    return seq[p - 1] if 1 <= p <= len(seq) else "N"



# ---------------------------------------------------------
# STEP 8: Add nucleotide identity at each breakpoint
# ---------------------------------------------------------
# This step checks which base (A/C/G/T) is present at every
# breakpoint. It also summarizes nucleotide counts by strand.
# This helps check if breakpoints are enriched near specific bases.
# ---------------------------------------------------------
def add_breakpoint_bases(output_folder: Path) -> tuple[Path, Path]:
    # add nucleotide at each breakpoint
    df = pd.read_csv(output_folder / "plus_minus_MFE.csv")
    plus_seq = read_fasta_sequence(output_folder / "payload_plusstrand.fasta")
    minus_seq = read_fasta_sequence(output_folder / "payload_minusstrand.fasta")

    df["Breakpoint"] = df["Breakpoint"].astype(int)
    df["Strand"] = df["Category"].apply(
        lambda x: "plus" if str(x).startswith("Coding") else "minus"
    )
    df["Base"] = df.apply(
        lambda row: get_base_1based(
            plus_seq if row["Strand"] == "plus" else minus_seq,
            row["Breakpoint"],
        ),
        axis=1,
    )

    annotated_path = output_folder / "plus_minus_MFE_with_bases.csv"
    df.to_csv(annotated_path, index=False)

    strand_counts = []
    valid_bases = {"A", "C", "G", "T"}
    for strand, subdf in df.groupby("Strand"):
        counts = subdf.groupby("Base", as_index=False)["Count"].sum()
        counts["Base"] = counts["Base"].str.upper()
        counts.loc[~counts["Base"].isin(valid_bases), "Base"] = "N"
        counts = counts.groupby("Base", as_index=False)["Count"].sum()
        total = counts["Count"].sum()
        counts["Proportion"] = counts["Count"] / total if total else 0
        counts.insert(0, "Strand", strand)
        strand_counts.append(counts)

    final_counts = pd.concat(strand_counts, ignore_index=True) if strand_counts else pd.DataFrame()
    counts_path = output_folder / "nucleotide_counts_by_strand.csv"
    final_counts.to_csv(counts_path, index=False)

    print(f"done: Wrote annotated breakpoint table: {annotated_path}")
    print(f"done: Wrote nucleotide count table:    {counts_path}")
    return annotated_path, counts_path



# ---------------------------------------------------------
# Helper: calculate payload background base frequency
# ---------------------------------------------------------
# This gives the expected A/C/G/T frequency from the full payload
# sequence, used as the background for enrichment testing.
# ---------------------------------------------------------
def background_freq(seq: str) -> dict[str, float]:
    total = len(seq)
    if total == 0:
        return {base: 0.0 for base in "ACGT"}
    return {base: seq.count(base) / total for base in "ACGT"}



# ---------------------------------------------------------
# STEP 9: Nucleotide enrichment test
# ---------------------------------------------------------
# This compares observed breakpoint base counts to expected base
# counts based on the payload background sequence. The chi-square
# test is used to check whether breakpoint bases differ from the
# overall payload nucleotide composition.
# ---------------------------------------------------------
def run_chi_square(output_folder: Path) -> tuple[Path, Path]:
    # compare breakpoint bases with payload background bases
    count_file = output_folder / "nucleotide_counts_by_strand.csv"
    df = pd.read_csv(count_file)
    df["Base"] = df["Base"].str.upper()
    df = df[df["Base"].isin(list("ACGT"))]

    plus_seq = read_fasta_sequence(output_folder / "payload_plusstrand.fasta")
    minus_seq = read_fasta_sequence(output_folder / "payload_minusstrand.fasta")
    bg_by_strand = {
        "plus": background_freq(plus_seq),
        "minus": background_freq(minus_seq),
    }

    rows_summary = []
    rows_detailed = []
    bases = ["A", "C", "G", "T"]

    for strand, sub in df.groupby("Strand"):
        observed = np.array([sub.loc[sub["Base"] == b, "Count"].sum() for b in bases], dtype=float)
        total_obs = observed.sum()
        bg = bg_by_strand[strand.lower()]
        expected = np.array([bg[b] * total_obs for b in bases], dtype=float)

        if total_obs == 0 or np.any(expected == 0):
            chi2, p_value = np.nan, np.nan
        else:
            chi2, p_value = chisquare(f_obs=observed, f_exp=expected)

        with np.errstate(divide="ignore", invalid="ignore"):
            std_resid = (observed - expected) / np.sqrt(expected)
            std_resid = np.where(expected > 0, std_resid, np.nan)

        rows_summary.append(
            {
                "Strand": strand,
                "TotalBreakpoints": int(total_obs),
                "Chi2": round(float(chi2), 3) if not np.isnan(chi2) else np.nan,
                "p_value": float(p_value) if not np.isnan(p_value) else np.nan,
                "bg_A": round(bg["A"], 4),
                "bg_C": round(bg["C"], 4),
                "bg_G": round(bg["G"], 4),
                "bg_T": round(bg["T"], 4),
            }
        )

        for i, base in enumerate(bases):
            rows_detailed.append(
                {
                    "Strand": strand,
                    "Base": base,
                    "Observed": int(observed[i]),
                    "Expected": round(float(expected[i]), 3),
                    "StdResidual": round(float(std_resid[i]), 3),
                }
            )

    summary_path = output_folder / "chi_square_results_by_strand.csv"
    detailed_path = output_folder / "chi_square_detailed_by_strand.csv"
    pd.DataFrame(rows_summary).to_csv(summary_path, index=False)
    pd.DataFrame(rows_detailed).to_csv(detailed_path, index=False)

    print(f"done: Wrote chi-square summary:  {summary_path}")
    print(f"done: Wrote chi-square details:  {detailed_path}")
    return summary_path, detailed_path



# ---------------------------------------------------------
# STEP 10: Plot breakpoint distribution
# ---------------------------------------------------------
# This makes count and proportion bar plots for plus and minus
# strand breakpoints. Bars are colored by MFE value so breakpoint
# hotspots can be compared with predicted local structure.
# ---------------------------------------------------------
def plot_breakpoints(output_folder: Path) -> tuple[Path, Path]:
    # simple plots for count and proportion
    df = pd.read_csv(output_folder / "plus_minus_MFE.csv")
    df["Breakpoint"] = df["Breakpoint"].astype(int)

    plot_files = []
    for y_col, ylabel, out_name in [
        ("Proportion", "Proportion", "Figure_Proportion.png"),
        ("Count", "Count", "Figure_Count.png"),
    ]:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
        categories = ["Coding / Plus strand", "Non-Coding / Minus strand"]

        for ax, category in zip(axes, categories):
            sub = df[df["Category"] == category].sort_values("Breakpoint")
            bars = ax.bar(sub["Breakpoint"], sub[y_col], width=20)

            if not sub.empty:
                mfe_values = sub["MFE"].to_numpy(dtype=float)
                norm = plt.Normalize(vmin=np.nanmin(mfe_values), vmax=np.nanmax(mfe_values))
                cmap = plt.cm.viridis
                for bar, mfe in zip(bars, mfe_values):
                    bar.set_color(cmap(norm(mfe)))

            ax.set_title(category)
            ax.set_xlabel("Snapback coordinate")
            ax.tick_params(axis="x", labelrotation=90, labelsize=8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        axes[0].set_ylabel(ylabel)
        fig.suptitle(f"Breakpoint Distribution ({ylabel})")
        fig.tight_layout()

        out_path = output_folder / out_name
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        fig.savefig(output_folder / out_name.replace(".png", ".pdf"), bbox_inches="tight")
        plt.close(fig)
        plot_files.append(out_path)
        print(f"done: Wrote plot: {out_path}")

    return tuple(plot_files)



# ---------------------------------------------------------
# FINAL PIPELINE: Run all steps for one sample
# ---------------------------------------------------------
# This function runs the complete workflow in order. Keeping it in
# one function makes it easy to loop over many samples without
# copying and pasting the same notebook cells again and again.
# ---------------------------------------------------------
def run_one_sample(tile_file: Path, reference_fasta: Path, output_folder: Path) -> None:
    # run everything for one sample
    tile_file = Path(tile_file)
    reference_fasta = Path(reference_fasta)
    output_folder = Path(output_folder)

    if not tile_file.exists():
        raise FileNotFoundError(f"Tile file not found: {tile_file}")
    if not reference_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")

    if RESET_OUTPUT_FOLDER and output_folder.exists():
        shutil.rmtree(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 80)
    print("Running sample")
    print(f"Tile file:       {tile_file}")
    print(f"Reference FASTA: {reference_fasta}")
    print(f"Output folder:   {output_folder}")
    print("=" * 80)

    plus_folder = output_folder / "plus"
    minus_folder = output_folder / "minus"
    plus_folder.mkdir(parents=True, exist_ok=True)
    minus_folder.mkdir(parents=True, exist_ok=True)

    parsed_csv = parse_snapback_tile_file(tile_file, output_folder)
    plus_csv, minus_csv = summarize_breakpoints(parsed_csv, output_folder)
    plus_fasta, minus_fasta = write_payload_fastas(reference_fasta, output_folder)

    plus_windows = write_breakpoint_windows(
        plus_fasta,
        plus_csv,
        "EndA",
        plus_folder / "Final_plus_Comb_brkp.seq",
    )
    minus_windows = write_breakpoint_windows(
        minus_fasta,
        minus_csv,
        "StartB",
        minus_folder / "Final_minus_Comb_brkp.seq",
    )

    plus_rnafold = run_rnafold(plus_windows, plus_folder / "Comb_plus_rnafold.seq")
    minus_rnafold = run_rnafold(minus_windows, minus_folder / "Comb_minus_rnafold.seq")

    extract_mfe_from_rnafold(
        plus_rnafold,
        "EndA",
        output_folder / "Drp_duplicates_column_of_Intrst_plus.csv",
    )
    extract_mfe_from_rnafold(
        minus_rnafold,
        "StartB",
        output_folder / "Drp_duplicates_column_of_Intrst_minus.csv",
    )

    merge_plus_minus_mfe(output_folder)
    add_breakpoint_bases(output_folder)
    run_chi_square(output_folder)
    plot_breakpoints(output_folder)

    print(f"done: Finished sample. Results are in: {output_folder}")



# ---------------------------------------------------------
# Command-line options
# ---------------------------------------------------------
# This lets the same script run from terminal using --sample.
# Multiple --sample entries can be provided for batch processing.
# ---------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Snapback RNAfold analysis"
    )
    parser.add_argument(
        "--sample",
        nargs=3,
        action="append",
        metavar=("TILE_FILE", "REFERENCE_FASTA", "OUTPUT_FOLDER"),
        help="One sample: tile-count file, reference FASTA, and output folder. Repeat for multiple samples.",
    )
    parser.add_argument(
        "--rnafold",
        default=RNAFOLD_CMD,
        help="RNAfold command/path. Default: RNAfold",
    )
    parser.add_argument(
        "--keep-output-folder",
        action="store_true",
        help="Do not delete the output folder before running.",
    )
    return parser.parse_args()



# ---------------------------------------------------------
# Main function
# ---------------------------------------------------------
# This checks whether samples were given from the command line.
# If not, it uses FILE_REFERENCE_PAIRS from the top of the script.
# Then it loops through each sample and runs the full pipeline.
# ---------------------------------------------------------
def main() -> None:
    global RNAFOLD_CMD, RESET_OUTPUT_FOLDER

    args = parse_args()
    RNAFOLD_CMD = args.rnafold
    RESET_OUTPUT_FOLDER = not args.keep_output_folder

    if args.sample:
        sample_pairs = [(Path(a), Path(b), Path(c)) for a, b, c in args.sample]
    else:
        sample_pairs = FILE_REFERENCE_PAIRS

    if not sample_pairs:
        raise SystemExit(
            "No samples provided. Use --sample TILE_FILE REFERENCE_FASTA OUTPUT_FOLDER, "
            "or edit FILE_REFERENCE_PAIRS inside the script."
        )

    for tile_file, reference_fasta, output_folder in sample_pairs:
        run_one_sample(tile_file, reference_fasta, output_folder)


if __name__ == "__main__":
    main()
