# ATAC-seq workflow and current progress

## Step 0 — define comparison groups (recorded done)

- Young: PDL20 (A/B/C; three bigWig replicates).
- Late/RS: PDL46 + PDL50 (A/B; four bigWig replicates total).

Keep these groups fixed across all scripts unless a documented project decision changes them.

## Step 1 — signal over atlas peaks (recorded done)

Script: `2_scripts/1_Preprocessing/1_atac_atlas_with_signal.py`

The atlas is a union set of genomic intervals; a peak is not intrinsically “a PDL20 peak” or “a PDL50 peak.” The same coordinate is measured in each group. The output records group-level signal and the contrast between them.

Recorded input atlas: 363,470 peaks.

Recorded output:

`1_data_raw/1_preprocessing_and_intermediate/atac_atlas_with_signal.csv`

Recorded derived values: `young_mean`, `rs_mean`, and `log2FC_rs_vs_young` (with the other columns defined by the script).

## Step 2 — select 10,000 benchmark intervals (recorded done)

Script: `2_scripts/1_Preprocessing/2_Select_benchmark_regions.py`

Recorded selection:

- Filter to atlas peaks with `max_signal` at or above the 50th percentile (median; recorded cutoff 0.12655201).
- Select 2,500 RS-gained intervals by largest log2 fold change.
- Select 2,500 RS-lost intervals by smallest log2 fold change.
- Select 2,500 stable, accessible intervals with `abs(log2FC) <= 0.25`.
- Create 2,500 negative intervals by shuffling stable intervals outside the atlas, constrained by chromosome sizes and exclusions.

Output: `1_data_raw/2_benchmark/benchmark_10000.bed` with columns `chr`, `start`, `end`, `name`, and `label`.

**Caveat:** the shuffled negatives are constructed outside the atlas and may be easier to distinguish than the biological classes. Keep them as a separately reported control and avoid over-interpreting their performance.

## Step 3 — experimental ground truth for benchmark (recorded done)

Script: `2_scripts/1_Preprocessing/3_Benchmark_groundtruth_ATAC.py`

For each of the 10,000 regions, the script extracts mean bigWig signal across the three PDL20 replicates and the four PDL46/PDL50 replicates, then computes the logged contrast using the implementation's epsilon.

Recorded output: `1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv` (10,000 rows, 8 columns) and a completion marker.

Before reusing it, inspect the script for the exact output column definitions, epsilon, coordinate handling, missing-value behavior, and bigWig signal semantics. Record these in the results report.

## Step 4.0 — reference (recorded done)

The local UCSC hg38 FASTA and `.fai` index have been created under `1_data_raw/0_reference/hg38/`.

## Step 4A — AlphaGenome ATAC pilot (verify completion)

Script: `2_scripts/1_Preprocessing/4A_AlphaGenome_ATAC_pilot.py`

Expected outputs described in the script notes include an interval JSON, extracted sequence FASTA, ATAC track metadata, pilot scores, candidate fibroblast-track table, log, and completion marker. Confirm the actual script's output names and inspect its log; a log file in the directory alone does not prove the pilot succeeded.

Choose a track only after checking the metadata. Record whether it is actually WI-38, another fibroblast, or a more distant proxy. Do not call a proxy a matched WI-38 track.

## Step 4B — proposed batch benchmark

Before making 10,000 API calls:

1. Validate one region end-to-end and confirm interval coordinates, sequence length, output array shapes, requested track and signal aggregation.
2. Write down the evaluation metric and baseline before running all regions.
3. Implement rate-limit handling, retries, checkpointing, safe resume, and a clear completed/failed row status.
4. Save output locally under `1_data_raw/3_alphagenome/` and keep it out of Git.
5. Evaluate by `RS_gained`, `RS_lost`, `Stable`, and `Negative`, and inspect failure cases.

Avoid implying that the model predicts the young-to-senescent fold change unless the model is actually provided appropriate condition-specific input and that capability is validated. A single sequence-conditioned prediction can be compared with each experimental state, but it does not by itself observe the state label.
