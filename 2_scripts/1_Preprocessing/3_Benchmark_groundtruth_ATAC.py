"""
Step 3: Compute experimental ATAC ground truth for benchmark regions (10,000)

Reads benchmark_10000.be (chr start end name label) and extracts mean bigWig
signal per region across replicate samples:

Young group : PDL20 (3 replicates)
RS/Late group : PDL46 + PDL50 (4 replicates)

Then averages within group and computes:
log2FC = log2((rs_mean + eps)/(young_mean + eps))

Outputs:
- benchmarks_10000_groundtruth.csv
- benchmarks_10000_groundtruth.done.txt

"""

#-----------------------------Imports--------------------------------------------#
import os
import pandas as pd
import numpy as np
import pyBigWig
from tqdm import tqdm

from utils_log import setup_logger, write_done_marker
#--------------------------------------------------------------------------------#

#------------------------Data-Loading------------------------#
IN_DIR   = "/home/ibab/1_POJECTS/Telos_NGS_Project/1_data_raw"
OUT_DIR  = os.path.join(IN_DIR, "2_benchmark")
RAW_DIR  = os.path.join(IN_DIR, "GSE175533_RAW")

BENCH_BED = os.path.join(OUT_DIR, "benchmark_10000.bed")
OUT_CSV   = os.path.join(OUT_DIR, "benchmark_10000_groundtruth.csv")

DONE_MARKER = os.path.join(OUT_DIR, "benchmark_10000_groundtruth.done.txt")
LOG_DIR     = os.path.join(OUT_DIR, "logs")

EPS = 1e-3

young_filenames = [
    "GSM5341071_PDL20_TP1_omni_A_S4_R1_001.filter.rmdup.bw",
    "GSM5341072_PDL20_TP1_omni_B_S5_R1_001.filter.rmdup.bw",
    "GSM5341073_PDL20_TP1_omni_C_S6_R1_001.filter.rmdup.bw",
]
rs_filenames = [
    "GSM5341086_PDL46_TP6_omni_A_S42_R1_001.filter.rmdup.bw",
    "GSM5341087_PDL46_TP6_omni_B_S43_R1_001.filter.rmdup.bw",
    "GSM5341088_PDL50_TP7_omni_A_S48_R1_001.filter.rmdup.bw",
    "GSM5341089_PDL50_TP7_omni_B_S49_R1_001.filter.rmdup.bw",
]
#-----------------------------------------------------------#

def mean_signal(bw_path, regions):
    bw = pyBigWig.open(bw_path)
    out = np.empty(len(regions), dtype=np.float32)
    for i, (c, s, e) in enumerate(tqdm(regions, desc=os.path.basename(bw_path))):
        v = bw.stats(c, int(s), int(e), type="mean")[0]
        out[i] = 0.0 if (v is None or np.isnan(v)) else float(v)
    bw.close()
    return out

def main():
    logger, log_path = setup_logger("step3_benchmark_groundtruth_atac", LOG_DIR)

    # sanity checks
    if not os.path.exists(BENCH_BED):
        raise FileNotFoundError(BENCH_BED)

    young_bw = [os.path.join(RAW_DIR, f) for f in young_filenames]
    rs_bw    = [os.path.join(RAW_DIR, f) for f in rs_filenames]
    for p in young_bw + rs_bw:
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    logger.info(f"Loading benchmark bed: {BENCH_BED}")
    bench = pd.read_csv(
        BENCH_BED, sep="\t", header=None,
        names=["chr", "start", "end", "name", "label"]
    )
    logger.info(f"Benchmark rows: {bench.shape[0]}")

    regions = list(bench[["chr", "start", "end"]].itertuples(index=False, name=None))

    logger.info("Extracting signals: Young (PDL20) replicates")
    Y = np.mean([mean_signal(p, regions) for p in young_bw], axis=0)

    logger.info("Extracting signals: RS/Late (PDL46+PDL50) replicates")
    R = np.mean([mean_signal(p, regions) for p in rs_bw], axis=0)

    out = bench.copy()
    out["young_mean"] = Y
    out["rs_mean"]    = R
    out["log2FC_rs_vs_young"] = np.log2((out["rs_mean"] + EPS) / (out["young_mean"] + EPS))

    out.to_csv(OUT_CSV, index=False)
    logger.info(f"Wrote: {OUT_CSV} shape={out.shape}")

    done_text = (
        f"SUCCESS\n"
        f"OUT_CSV={OUT_CSV}\n"
        f"BENCH_BED={BENCH_BED}\n"
        f"YOUNG_N={len(young_bw)}\n"
        f"RS_N={len(rs_bw)}\n"
        f"LOG={log_path}\n"
    )
    write_done_marker(DONE_MARKER, done_text)
    logger.info(f"Wrote DONE marker: {DONE_MARKER}")
    logger.info("STEP 3 COMPLETE")


if __name__ == "__main__":
    main()
