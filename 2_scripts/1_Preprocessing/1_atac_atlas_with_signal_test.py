
"""
This script quantifies and compares ATAC-seq chromatin accessibility between
early passage (young, PDL20) and late passage (replicative senescent, PDL46–50)
cells.

It loads peak regions from atac_atlas.bed and uses pyBigWig to extract the mean
accessibility signal for each region across:
- 3 young replicate bigWigs (PDL20)
- 4 senescent/late replicate bigWigs (PDL46 + PDL50)

Then it averages signal across replicates and computes differential accessibility:
log2FC = log2( (RS + eps) / (Young + eps) ), eps=1e-3.

Output:
- atac_atlas_with_signal_test.csv
"""

#------------------------------------------#
import pyBigWig
import pandas as pd
import numpy as np
from tqdm import tqdm
import os
#------------------------------------------#


#------------------------Data-Loading------------------------#

# Resolve the project root relative to this script
PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)

IN_DIR = os.path.join(PROJECT_ROOT, "1_data_raw")

# Input BED file
BED_FILE = os.path.join(IN_DIR, "atac_atlas.bed")

bed = pd.read_csv(
    BED_FILE,
    sep="\t",
    header=None,
    names=["chr", "start", "end", "id"]
)

# IMPORTANT: make it a LIST, not an iterator
regions = list(
    bed[["chr", "start", "end"]].itertuples(
        index=False, name=None
    )
)

# Directory containing replicate bigWig files
data_dir = os.path.join(IN_DIR, "GSE175533_RAW")

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

young_bw = [
    os.path.join(data_dir, f) for f in young_filenames
]

rs_bw = [
    os.path.join(data_dir, f) for f in rs_filenames
]

# Output directory and file
OUTPUT_DIR = os.path.join(
    IN_DIR,
    "1_preprocessing_and_intermediate"
)

os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_FILE = os.path.join(
    OUTPUT_DIR,
    "atac_atlas_with_signal_test.csv"
)

#-----------------------------------------------------------#


#------------------------Signal Extraction------------------#

def mean_signal(bw_path, regions):

    bw = pyBigWig.open(bw_path)

    out = np.empty(len(regions), dtype=np.float32)

    for i, (c, s, e) in enumerate(
        tqdm(regions, desc=os.path.basename(bw_path))
    ):

        v = bw.stats(c, int(s), int(e), type="mean")[0]

        out[i] = (
            0.0 if (v is None or np.isnan(v)) else float(v)
        )

    bw.close()

    return out

#-----------------------------------------------------------#


#------------------------Main Analysis---------------------#

def main():

    Y = np.mean(
        [mean_signal(p, regions) for p in young_bw],
        axis=0
    )

    R = np.mean(
        [mean_signal(p, regions) for p in rs_bw],
        axis=0
    )

    out = bed.copy()

    out["young_mean"] = Y
    out["rs_mean"] = R

    eps = 1e-3

    out["log2FC_rs_vs_young"] = np.log2(
        (out["rs_mean"] + eps) /
        (out["young_mean"] + eps)
    )

    out.to_csv(OUTPUT_FILE, index=False)

    print("Wrote:", OUTPUT_FILE)
    print("Output shape:", out.shape)

#-----------------------------------------------------------#


if __name__ == "__main__":
    main()
