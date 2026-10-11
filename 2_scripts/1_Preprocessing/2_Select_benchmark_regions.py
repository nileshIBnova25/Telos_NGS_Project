"""
Step 2: Create a fixed benchmark region set (10,000 regions) from atlas peaks + signals.

Inputs:
- atac_atlas_with_signal.csv  (Step 1 output)
- GSE175533_atlas.bed.gz      (for exclusion when generating negatives)
- one bigWig file             (only used to export chromosome sizes)

Method:
- Filter to WI-38-relevant peaks using max_signal quantile filter.
- Select:
  - N RS_gained peaks (largest log2FC)
  - N RS_lost peaks   (smallest log2FC)
  - N Stable peaks    (abs(log2FC) <= stable_abs_log2fc among accessible peaks)
  - N Negatives       (bedtools shuffle stable intervals outside atlas)

Output:
- benchmark_10000.bed with columns: chr start end name label
- logs + .done marker
"""

#------------------------Imports-----------------------------#
import os
import gzip
import subprocess
import pandas as pd
import numpy as np
import pyBigWig

from utils_log import setup_logger, write_done_marker
#-----------------------------------------------------------#


#------------------------Data-Loading------------------------#
IN_DIR   = "/home/ibab/1_POJECTS/Telos_NGS_Project/1_data_raw"
OUT_DIR  = os.path.join(IN_DIR, "2_benchmark")
RAW_DIR  = os.path.join(IN_DIR, "GSE175533_RAW")

IN_CSV   = os.path.join(IN_DIR, "1_preprocessing_and_intermediate/atac_atlas_with_signal.csv")
ATLAS_GZ = os.path.join(IN_DIR, "GSE175533_atlas.bed.gz")

CHROM_SIZES = os.path.join(OUT_DIR, "hg38.chrom.sizes.from_bw.txt")
ATLAS_EXCL  = os.path.join(OUT_DIR, "atlas_exclude.bed")

BW_FOR_CHROMS = os.path.join(RAW_DIR, "GSM5341071_PDL20_TP1_omni_A_S4_R1_001.filter.rmdup.bw")

N = 2500  # 4 bins * 2500 = 10000

# selection params (keep simple)
MAX_SIGNAL_KEEP_QUANTILE = 0.50   # keep top 50% peaks by max(young_mean, rs_mean)
STABLE_ABS_LOG2FC = 0.25

OUT_BED = os.path.join(OUT_DIR, "benchmark_10000.bed")
DONE_MARKER = os.path.join(OUT_DIR, "benchmark_10000.done.txt")
LOG_DIR = os.path.join(OUT_DIR, "logs")
#-----------------------------------------------------------#


def require_bedtools(logger):
    try:
        r = subprocess.run(["bedtools", "--version"], check=True, capture_output=True, text=True)
        logger.info("bedtools OK: " + r.stdout.strip())
    except Exception:
        raise RuntimeError("bedtools not found. Install with: conda install -c bioconda bedtools -y")


def make_chrom_sizes(logger):
    if os.path.exists(CHROM_SIZES):
        logger.info(f"Chrom sizes exists: {CHROM_SIZES}")
        return
    logger.info(f"Writing chrom sizes from: {BW_FOR_CHROMS}")
    bw = pyBigWig.open(BW_FOR_CHROMS)
    with open(CHROM_SIZES, "w") as f:
        for c, l in bw.chroms().items():
            f.write(f"{c}\t{l}\n")
    bw.close()
    logger.info(f"Wrote: {CHROM_SIZES}")


def make_atlas_exclude(logger):
    if os.path.exists(ATLAS_EXCL):
        logger.info(f"Atlas exclude exists: {ATLAS_EXCL}")
        return
    logger.info(f"Writing atlas exclude bed from: {ATLAS_GZ}")
    with gzip.open(ATLAS_GZ, "rt") as fin, open(ATLAS_EXCL, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            fout.write("\t".join(parts[:3]) + "\n")
    logger.info(f"Wrote: {ATLAS_EXCL}")


def write_bed(df, path):
    df[["chr","start","end","name","label"]].to_csv(path, sep="\t", header=False, index=False)


def main():
    logger, log_path = setup_logger("step2_make_benchmark_10000", LOG_DIR)
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---- sanity checks
    for p in [IN_CSV, ATLAS_GZ, BW_FOR_CHROMS]:
        if not os.path.exists(p):
            raise FileNotFoundError(p)
    require_bedtools(logger)

    # ---- prepare auxiliary files
    make_chrom_sizes(logger)
    make_atlas_exclude(logger)

    # ---- load signals
    logger.info(f"Loading CSV: {IN_CSV}")
    df = pd.read_csv(IN_CSV)
    logger.info(f"CSV rows: {df.shape[0]} cols: {df.shape[1]}")
    required_cols = {"chr","start","end","id","young_mean","rs_mean","log2FC_rs_vs_young"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in CSV: {missing}")

    # ---- filter to WI-38-relevant peaks
    df["max_signal"] = df[["young_mean","rs_mean"]].max(axis=1)
    cut = df["max_signal"].quantile(MAX_SIGNAL_KEEP_QUANTILE)
    df2 = df[df["max_signal"] >= cut].copy()
    logger.info(f"Filter max_signal >= q{MAX_SIGNAL_KEEP_QUANTILE}: keep {len(df2)} peaks (cut={cut})")

    # ---- select gained/lost
    gained = df2.sort_values("log2FC_rs_vs_young", ascending=False).head(N).copy()
    gained["label"] = "RS_gained"
    gained["name"] = gained["label"] + "_" + gained["id"].astype(str)

    lost = df2.sort_values("log2FC_rs_vs_young", ascending=True).head(N).copy()
    lost["label"] = "RS_lost"
    lost["name"] = lost["label"] + "_" + lost["id"].astype(str)

    # ---- select stable
    stable_pool = df2[df2["log2FC_rs_vs_young"].abs() <= STABLE_ABS_LOG2FC].copy()
    logger.info(f"Stable pool size (abs(log2FC)<={STABLE_ABS_LOG2FC}): {len(stable_pool)}")
    if len(stable_pool) < N:
        raise RuntimeError(f"Not enough stable peaks ({len(stable_pool)}) for N={N}. "
                           f"Relax STABLE_ABS_LOG2FC or reduce N.")
    stable = stable_pool.sample(n=N, random_state=1).copy()
    stable["label"] = "Stable"
    stable["name"] = stable["label"] + "_" + stable["id"].astype(str)

    # ---- write intermediate beds
    gained_bed = os.path.join(OUT_DIR, "bench_gained.bed")
    lost_bed   = os.path.join(OUT_DIR, "bench_lost.bed")
    stable_bed = os.path.join(OUT_DIR, "bench_stable.bed")
    write_bed(gained, gained_bed)
    write_bed(lost, lost_bed)
    write_bed(stable, stable_bed)
    logger.info(f"Wrote: {gained_bed}, {lost_bed}, {stable_bed}")

    # ---- negatives via bedtools shuffle outside atlas
    stable_seed = os.path.join(OUT_DIR, "stable_seed_for_negatives.bed")
    stable[["chr","start","end"]].to_csv(stable_seed, sep="\t", header=False, index=False)

    neg_tmp = os.path.join(OUT_DIR, "negatives_tmp.bed")
    cmd = [
        "bedtools", "shuffle",
        "-i", stable_seed,
        "-g", CHROM_SIZES,
        "-excl", ATLAS_EXCL,
        "-chrom",
        "-maxTries", "10000",
    ]
    logger.info("Running: " + " ".join(cmd))
    with open(neg_tmp, "w") as f:
        subprocess.run(cmd, check=True, stdout=f)

    neg = pd.read_csv(neg_tmp, sep="\t", header=None, names=["chr","start","end"])
    if len(neg) < N:
        raise RuntimeError(f"bedtools shuffle produced {len(neg)} negatives; need {N}. Increase maxTries.")
    neg = neg.head(N).copy()
    neg["name"]  = ["NEG_" + str(i+1) for i in range(len(neg))]
    neg["label"] = "Negative"

    neg_bed = os.path.join(OUT_DIR, "bench_negative.bed")
    write_bed(neg, neg_bed)
    logger.info(f"Wrote: {neg_bed}")

    # ---- combine final benchmark bed
    all_df = pd.concat([
        gained[["chr","start","end","name","label"]],
        lost[["chr","start","end","name","label"]],
        stable[["chr","start","end","name","label"]],
        neg[["chr","start","end","name","label"]],
    ], axis=0)

    all_df.to_csv(OUT_BED, sep="\t", header=False, index=False)
    logger.info(f"Wrote final benchmark bed: {OUT_BED}")
    logger.info("Counts by label:\n" + str(all_df["label"].value_counts()))

    # ---- done marker for collaborators
    done_text = (
        f"SUCCESS\n"
        f"OUT_BED={OUT_BED}\n"
        f"N_PER_CLASS={N}\n"
        f"TOTAL={len(all_df)}\n"
        f"IN_CSV={IN_CSV}\n"
        f"ATLAS_GZ={ATLAS_GZ}\n"
        f"LOG={log_path}\n"
    )
    write_done_marker(DONE_MARKER, done_text)
    logger.info(f"Wrote DONE marker: {DONE_MARKER}")
    logger.info("STEP 2 COMPLETE")


if __name__ == "__main__":
    main()
