"""
step 4A: AlphaGenome ATAC pilot + track inventory (with local hg38 sequence saved)

Input:
- benchmark_10000.bed
- hg38.chrom.sizes.from_bw.txt
- local hg38 FASTA (UCSC chr naming)
- AlphaGenome API key from env : ALPHA_GENOME_API_KEY or ALPHAGENOME_API_KEY

Output (in 1_data_raw/3_alphagenome/):
- pilot_interval.json
- pilot_sequence.fa
- pilot_atac_track_metadata.csv
- pilot_atac_track_scores.csv
- pilot_fibroblast_candidate_tracks.csv
- step4A_alphagenome_atac_pilot.done.txt
- logs/step4A_*.log
"""

#----------------------Imports-----------------------------#
import os, json
import pandas as pd
import numpy as np
from pyfaidx import Fasta

from utils_log import setup_logger, write_done_marker

from alphagenome.models import dna_client
from alphagenome.data import genome
#----------------------------------------------------------#

#------------------------Data-Loading------------------------#
IN_DIR = "/home/ibab/1_POJECTS/Telos_NGS_Project/1_data_raw"

BENCHMARK_DIR = os.path.join(IN_DIR, "2_benchmark")

AG_BASE_DIR = os.path.join(IN_DIR, "3_alphagenome")
OUT_DIR = os.path.join(AG_BASE_DIR, "00_pilot")
LOG_DIR = os.path.join(OUT_DIR, "logs")

BENCH_BED = os.path.join(BENCHMARK_DIR, "benchmark_10000.bed")
CHROM_SIZES = os.path.join(BENCHMARK_DIR, "hg38.chrom.sizes.from_bw.txt")

HG38_FASTA = os.path.join(IN_DIR, "0_reference/hg38/hg38.fa")
#-----------------------------------------------------------#
PILOT_LABEL = "RS_gained"
SEQUENCE_LENGTH = 16384  # pilot context

METADATA_CSV   = os.path.join(OUT_DIR, "pilot_atac_track_metadata.csv")
SCORES_CSV     = os.path.join(OUT_DIR, "pilot_atac_track_scores.csv")
CANDIDATES_CSV = os.path.join(OUT_DIR, "pilot_fibroblast_candidate_tracks.csv")
INTERVAL_JSON  = os.path.join(OUT_DIR, "pilot_interval.json")
PILOT_FASTA    = os.path.join(OUT_DIR, "pilot_sequence.fa")

DONE_MARKER = os.path.join(OUT_DIR, "step4A_alphagenome_atac_pilot.done.txt")
#-----------------------------------------------------------#

def get_api_key():
    key = os.getenv("ALPHA_GENOME_API_KEY") or os.getenv("ALPHAGENOME_API_KEY")
    if not key:
        raise RuntimeError(
            "AlphaGenome API key not found. Export ALPHA_GENOME_API_KEY (or ALPHAGENOME_API_KEY)."
        )
    return key

def load_chrom_sizes(path):
    chrom_sizes = {}
    with open(path) as f:
        for line in f:
            c, l = line.rstrip("\n").split("\t")[:2]
            chrom_sizes[c] = int(l)
    return chrom_sizes

def make_centered_interval(chrom, peak_start, peak_end, chrom_length, seq_len):
    midpoint = (int(peak_start) + int(peak_end)) // 2
    start = midpoint - (seq_len // 2)
    end = start + seq_len

    if start < 0:
        start = 0
        end = seq_len

    if end > chrom_length:
        end = chrom_length
        start = end - seq_len

    if start < 0 or (end - start) != seq_len:
        raise ValueError(f"Cannot build {seq_len}bp interval for {chrom}:{peak_start}-{peak_end}")

    return int(start), int(end)

def find_atac_output_type():
    OT = dna_client.OutputType
    for name in ["ATAC", "ATAC_SEQ", "ATACSEQ"]:
        if hasattr(OT, name):
            return getattr(OT, name)
    avail = [x for x in dir(OT) if x.isupper()]
    raise RuntimeError(f"Could not find ATAC OutputType. Available: {avail}")

def metadata_to_df(metadata):
    if isinstance(metadata, pd.DataFrame):
        return metadata.reset_index(drop=True)
    if hasattr(metadata, "to_pandas"):
        return metadata.to_pandas().reset_index(drop=True)
    if isinstance(metadata, (list, tuple)):
        rows = []
        for item in metadata:
            rows.append(vars(item) if hasattr(item, "__dict__") else item)
        return pd.DataFrame(rows)
    return pd.DataFrame(metadata)

def orient_values(values, n_tracks):
    v = np.asarray(values)
    v = np.squeeze(v)

    if v.ndim != 2:
        raise ValueError(f"ATAC values should be 2D; got shape {v.shape}")

    # positions x tracks
    if v.shape[1] == n_tracks:
        return v
    # tracks x positions
    if v.shape[0] == n_tracks:
        return v.T

    raise ValueError(f"Cannot match values shape {v.shape} with n_tracks={n_tracks}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    logger, log_path = setup_logger("step4A_alphagenome_atac_pilot", LOG_DIR)

    # sanity checks
    for p in [BENCH_BED, CHROM_SIZES, HG38_FASTA]:
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    api_key = get_api_key()
    logger.info("API key detected: YES (will not be printed)")

    bench = pd.read_csv(
        BENCH_BED, sep="\t", header=None,
        names=["chr","start","end","name","label"]
    )
    pool = bench[bench["label"] == PILOT_LABEL]
    if pool.empty:
        raise RuntimeError(f"No region found with label={PILOT_LABEL}")
    pilot = pool.iloc[0]

    chrom_sizes = load_chrom_sizes(CHROM_SIZES)
    chrom = pilot["chr"]
    if chrom not in chrom_sizes:
        raise RuntimeError(f"{chrom} not in chrom sizes file")

    context_start, context_end = make_centered_interval(
        chrom, pilot["start"], pilot["end"], chrom_sizes[chrom], SEQUENCE_LENGTH
    )

    # local sequence save (for other models)
    ref = Fasta(HG38_FASTA, rebuild=False)
    seq = str(ref[chrom][context_start:context_end]).upper()
    with open(PILOT_FASTA, "w") as f:
        f.write(f">{chrom}:{context_start}-{context_end}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
    logger.info(f"Wrote pilot local sequence FASTA: {PILOT_FASTA}")

    interval = genome.Interval(
        chromosome=str(chrom),
        start=int(context_start),
        end=int(context_end)
    )
    logger.info(f"Pilot peak: {chrom}:{pilot['start']}-{pilot['end']} ({pilot['name']})")
    logger.info(f"AlphaGenome interval: {chrom}:{context_start}-{context_end} (len={SEQUENCE_LENGTH})")

    atac_output = find_atac_output_type()
    model = dna_client.create(api_key)

    # IMPORTANT: ontology_terms is required; pass None for pilot
    output = model.predict_interval(
        interval=interval,
        requested_outputs=[atac_output],
        ontology_terms=None
    )

    if not hasattr(output, "atac") or output.atac is None:
        raise RuntimeError("AlphaGenome response has no 'atac' field.")

    atac = output.atac
    meta = metadata_to_df(atac.metadata)
    vals = orient_values(atac.values, n_tracks=len(meta))

    logger.info(f"ATAC values matrix: {vals.shape} (positions x tracks)")
    logger.info(f"ATAC metadata rows: {meta.shape[0]}")
    logger.info("Metadata columns: " + ", ".join(map(str, meta.columns)))

    # score predicted ATAC over the original peak span
    n_pos = vals.shape[0]
    rel_start = int(pilot["start"]) - context_start
    rel_end   = int(pilot["end"]) - context_start

    first_bin = int(np.floor(rel_start * n_pos / SEQUENCE_LENGTH))
    last_bin  = int(np.ceil(rel_end   * n_pos / SEQUENCE_LENGTH))
    first_bin = max(0, first_bin)
    last_bin  = min(n_pos, max(first_bin + 1, last_bin))

    track_scores = np.nanmean(vals[first_bin:last_bin, :], axis=0)

    meta = meta.reset_index(drop=True)
    meta["track_index"] = np.arange(len(meta))
    meta["pilot_pred_atac_mean"] = track_scores

    meta.to_csv(METADATA_CSV, index=False)
    meta.to_csv(SCORES_CSV, index=False)
    logger.info(f"Wrote: {METADATA_CSV}")
    logger.info(f"Wrote: {SCORES_CSV}")

    # candidate track filtering (depends on metadata fields; this is text-based and robust)
    text = meta.astype(str).agg(" | ".join, axis=1)
    pattern = r"WI[- ]?38|fibroblast|IMR[- ]?90|HFF|lung"
    candidates = meta[text.str.contains(pattern, case=False, regex=True, na=False)].copy()
    candidates.to_csv(CANDIDATES_CSV, index=False)
    logger.info(f"Candidate fibroblast/lung tracks: {len(candidates)} -> {CANDIDATES_CSV}")

    manifest = {
        "benchmark_name": str(pilot["name"]),
        "benchmark_label": str(pilot["label"]),
        "peak": f"{chrom}:{int(pilot['start'])}-{int(pilot['end'])}",
        "interval": f"{chrom}:{context_start}-{context_end}",
        "sequence_length": int(SEQUENCE_LENGTH),
        "hg38_fasta": HG38_FASTA,
        "metadata_csv": METADATA_CSV,
        "candidates_csv": CANDIDATES_CSV,
    }
    with open(INTERVAL_JSON, "w") as f:
        json.dump(manifest, f, indent=2)
    logger.info(f"Wrote: {INTERVAL_JSON}")

    done_text = (
        "SUCCESS\n"
        f"PILOT_PEAK={manifest['peak']}\n"
        f"INPUT_INTERVAL={manifest['interval']}\n"
        f"SEQUENCE_LENGTH={SEQUENCE_LENGTH}\n"
        f"ATAC_TRACKS={len(meta)}\n"
        f"CANDIDATE_TRACKS={len(candidates)}\n"
        f"PILOT_FASTA={PILOT_FASTA}\n"
        f"METADATA={METADATA_CSV}\n"
        f"CANDIDATES={CANDIDATES_CSV}\n"
        f"LOG={log_path}\n"
    )
    write_done_marker(DONE_MARKER, done_text)
    logger.info(f"Wrote DONE marker: {DONE_MARKER}")
    logger.info("STEP 4A COMPLETE")

if __name__ == "__main__":
    main()