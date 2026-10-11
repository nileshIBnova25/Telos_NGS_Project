"""
Step 4B: AlphaGenome ATAC batch prediction for benchmark_10000.bed (10,000 regions)

Robust version:
- single-instance lock to prevent accidental double runs
- frequent flushing to disk
- retry with exponential backoff
- safe stop handling

Input:
- 1_data_raw/2_benchmark/benchmark_10000.bed
- 1_data_raw/2_benchmark/hg38.chrom.sizes.from_bw.txt
- 1_data_raw/3_alphagenome/01_config/selected_atac_track.json
- AlphaGenome API key in env: ALPHA_GENOME_API_KEY (or ALPHAGENOME_API_KEY)

Output:
- 1_data_raw/3_alphagenome/02_atac_run/benchmark_10000_alphagenome_atac_pred.csv
- 1_data_raw/3_alphagenome/02_atac_run/step4B_alphagenome_atac_batch.done.txt
- 1_data_raw/3_alphagenome/02_atac_run/logs/step4B_*.log
"""

#----------------------Imports-----------------------------#
import os, json, time, signal, atexit
import numpy as np
import pandas as pd

from pathlib import Path
import sys
sys.path.append(str(Path(__file__).resolve().parents[1] / "1_Preprocessing"))
from utils_log import setup_logger, write_done_marker

from alphagenome.models import dna_client
from alphagenome.data import genome

import fcntl
#----------------------------------------------------------#

#------------------------Data-Loading------------------------#
IN_DIR = "/home/ibab/1_POJECTS/Telos_NGS_Project/1_data_raw"
BENCHMARK_DIR = os.path.join(IN_DIR, "2_benchmark")

AG_BASE_DIR = os.path.join(IN_DIR, "3_alphagenome")
CONFIG_DIR  = os.path.join(AG_BASE_DIR, "01_config")
OUT_DIR     = os.path.join(AG_BASE_DIR, "02_atac_run")
LOG_DIR     = os.path.join(OUT_DIR, "logs")

BENCH_BED   = os.path.join(BENCHMARK_DIR, "benchmark_10000.bed")
CHROM_SIZES = os.path.join(BENCHMARK_DIR, "hg38.chrom.sizes.from_bw.txt")
TRACK_JSON  = os.path.join(CONFIG_DIR, "selected_atac_track.json")

OUT_CSV     = os.path.join(OUT_DIR, "benchmark_10000_alphagenome_atac_pred.csv")
DONE_MARKER = os.path.join(OUT_DIR, "step4B_alphagenome_atac_batch.done.txt")

LOCK_PATH   = os.path.join(OUT_DIR, "step4B.lock")   # prevents 2 runs at once

SEQUENCE_LENGTH = 16384
MAX_RETRIES = 6
FLUSH_EVERY = 20          # flush often so you don't lose progress
SLEEP_BETWEEN_CALLS = 0.0 # set to 0.1 if you see throttling
#-----------------------------------------------------------#

STOP_REQUESTED = False
LOCK_FH = None

def _handle_stop(signum, frame):
    global STOP_REQUESTED
    STOP_REQUESTED = True

def get_api_key():
    key = os.getenv("ALPHA_GENOME_API_KEY") or os.getenv("ALPHAGENOME_API_KEY")
    if not key:
        raise RuntimeError("AlphaGenome API key not found in environment.")
    return key

def load_chrom_sizes(path):
    d = {}
    with open(path) as f:
        for line in f:
            c, l = line.rstrip("\n").split("\t")[:2]
            d[c] = int(l)
    return d

def make_centered_interval(chrom, start, end, chrom_len, seq_len):
    mid = (int(start) + int(end)) // 2
    s = mid - (seq_len // 2)
    e = s + seq_len
    if s < 0:
        s = 0
        e = seq_len
    if e > chrom_len:
        e = chrom_len
        s = e - seq_len
    if s < 0 or (e - s) != seq_len:
        raise ValueError(f"Bad interval for {chrom}:{start}-{end}")
    return int(s), int(e)

def find_atac_output_type():
    OT = dna_client.OutputType
    for name in ["ATAC", "ATAC_SEQ", "ATACSEQ"]:
        if hasattr(OT, name):
            return getattr(OT, name)
    avail = [x for x in dir(OT) if x.isupper()]
    raise RuntimeError(f"ATAC OutputType not found. Available: {avail}")

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
        raise ValueError(f"ATAC values must be 2D; got {v.shape}")
    if v.shape[1] == n_tracks:
        return v
    if v.shape[0] == n_tracks:
        return v.T
    raise ValueError(f"Cannot orient values shape={v.shape} vs n_tracks={n_tracks}")

def load_selected_track():
    with open(TRACK_JSON) as f:
        obj = json.load(f)
    if "track_index" not in obj:
        raise ValueError("selected_atac_track.json must contain 'track_index'")
    return int(obj["track_index"]), obj.get("note", "")

def acquire_lock(logger):
    """
    Prevent accidental second run writing to the same OUT_CSV.
    Uses an OS file lock. If already locked, exits immediately.
    """
    global LOCK_FH
    LOCK_FH = open(LOCK_PATH, "w")
    try:
        fcntl.flock(LOCK_FH, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError(f"Another Step4B run is already active (lock busy): {LOCK_PATH}")

    LOCK_FH.write(f"pid={os.getpid()}\n")
    LOCK_FH.flush()
    logger.info(f"Acquired lock: {LOCK_PATH}")

def release_lock():
    global LOCK_FH
    try:
        if LOCK_FH:
            fcntl.flock(LOCK_FH, fcntl.LOCK_UN)
            LOCK_FH.close()
            LOCK_FH = None
    except Exception:
        pass

def flush_rows(logger, rows):
    if not rows:
        return
    df_out = pd.DataFrame(rows)
    write_header = not os.path.exists(OUT_CSV)
    df_out.to_csv(OUT_CSV, mode="a", header=write_header, index=False)
    logger.info(f"Flushed {len(rows)} rows -> {OUT_CSV}")

def main():
    global STOP_REQUESTED

    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(CONFIG_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    logger, log_path = setup_logger("step4B_alphagenome_atac_batch", LOG_DIR)

    # stop handlers
    signal.signal(signal.SIGINT, _handle_stop)
    signal.signal(signal.SIGTERM, _handle_stop)
    atexit.register(release_lock)

    # sanity checks
    for p in [BENCH_BED, CHROM_SIZES, TRACK_JSON]:
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    acquire_lock(logger)

    api_key = get_api_key()
    track_index, track_note = load_selected_track()
    logger.info(f"Selected track_index={track_index} note={track_note}")

    chrom_sizes = load_chrom_sizes(CHROM_SIZES)
    atac_output = find_atac_output_type()
    model = dna_client.create(api_key)

    bench = pd.read_csv(BENCH_BED, sep="\t", header=None,
                        names=["chr", "start", "end", "name", "label"])

    # resume support
    done_names = set()
    if os.path.exists(OUT_CSV):
        prev = pd.read_csv(OUT_CSV)
        if "name" in prev.columns:
            done_names = set(prev["name"].astype(str))
        logger.info(f"Resume: existing rows={len(done_names)}")

    out_rows = []
    failed = 0
    cached_n_tracks = None

    for i, row in enumerate(bench.itertuples(index=False), start=1):
        if STOP_REQUESTED:
            logger.warning("STOP requested. Flushing remaining rows and exiting.")
            break

        name = str(row.name)
        if name in done_names:
            continue

        chrom = str(row.chr)
        if chrom not in chrom_sizes:
            logger.warning(f"Skipping {name}: chromosome {chrom} not in chrom sizes")
            continue

        ctx_s, ctx_e = make_centered_interval(chrom, row.start, row.end, chrom_sizes[chrom], SEQUENCE_LENGTH)
        interval = genome.Interval(chromosome=chrom, start=ctx_s, end=ctx_e)

        ok = False
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                output = model.predict_interval(
                    interval=interval,
                    requested_outputs=[atac_output],
                    ontology_terms=None
                )

                atac = output.atac
                meta = metadata_to_df(atac.metadata)
                vals = orient_values(atac.values, n_tracks=len(meta))

                if cached_n_tracks is None:
                    cached_n_tracks = vals.shape[1]
                    logger.info(f"Detected n_tracks={cached_n_tracks} on first success")

                if track_index < 0 or track_index >= vals.shape[1]:
                    raise IndexError(f"track_index {track_index} out of bounds for {vals.shape[1]} tracks")

                rel_s = int(row.start) - ctx_s
                rel_e = int(row.end) - ctx_s
                rel_s = max(0, rel_s)
                rel_e = min(SEQUENCE_LENGTH, rel_e)
                if rel_e <= rel_s:
                    rel_s, rel_e = 0, SEQUENCE_LENGTH

                pred_mean = float(np.nanmean(vals[rel_s:rel_e, track_index]))

                out_rows.append({
                    "name": name,
                    "label": str(row.label),
                    "chr": chrom,
                    "start": int(row.start),
                    "end": int(row.end),
                    "track_index": int(track_index),
                    "pred_atac_mean": pred_mean
                })
                ok = True
                break

            except Exception as e:
                wait = min(60.0, (2.0 ** attempt))  # exponential backoff up to 60s
                logger.warning(f"{name} attempt {attempt}/{MAX_RETRIES} failed: {repr(e)}; sleeping {wait:.1f}s")
                time.sleep(wait)

        if not ok:
            failed += 1

        if len(out_rows) >= FLUSH_EVERY:
            flush_rows(logger, out_rows)
            out_rows = []

        if SLEEP_BETWEEN_CALLS > 0:
            time.sleep(SLEEP_BETWEEN_CALLS)

    # final flush
    flush_rows(logger, out_rows)

    # done marker only if we reached end without STOP_REQUESTED
    if not STOP_REQUESTED:
        done_text = (
            "SUCCESS\n"
            f"OUT_CSV={OUT_CSV}\n"
            f"TRACK_INDEX={track_index}\n"
            f"FAILED={failed}\n"
            f"LOG={log_path}\n"
        )
        write_done_marker(DONE_MARKER, done_text)
        logger.info(f"Wrote DONE marker: {DONE_MARKER}")
        logger.info("STEP 4B COMPLETE")
    else:
        logger.warning("STOPPED before completion (no DONE marker written). You can rerun to resume.")

if __name__ == "__main__":
    main()