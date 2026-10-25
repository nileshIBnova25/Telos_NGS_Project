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
- atac_atlas_with_signal.csv
"""
#------------------------------------------#
import pyBigWig, pandas as pd, numpy as np
from tqdm import tqdm
import os
#------------------------------------------#
#--------------------------------------------------------------------------------#
# Bed file loading
bed = pd.read_csv("/home/ibab/1_POJECTS/Telos_NGS_Project/1_data_raw/atac_atlas.bed",sep="\t",header=None,names=["chr","start","end","id"])

# IMPORTANT: make it a LIST, not an iterator
regions = list(bed[["chr","start","end"]].itertuples(index=False, name=None))

data_dir = "/home/ibab/1_POJECTS/Telos_NGS_Project/1_data_raw/GSE175533_RAW"

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
young_bw = [os.path.join(data_dir,f) for f in young_filenames]
rs_bw = [os.path.join(data_dir,f) for f in rs_filenames]

#--------------------------------------------------------------------------------#
def mean_signal(bw_path, regions):
    bw = pyBigWig.open(bw_path)
    out = np.empty(len(regions), dtype=np.float32)
    for i, (c,s,e) in enumerate(tqdm(regions, desc=bw_path.split("/")[-1])):
        v = bw.stats(c, int(s), int(e), type="mean")[0]
        out[i] = 0.0 if (v is None or np.isnan(v)) else float(v)
    bw.close()
    return out
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
def main():
    Y = np.mean([mean_signal(p, regions) for p in young_bw], axis=0)
    R = np.mean([mean_signal(p, regions) for p in rs_bw], axis=0)

    out = bed.copy()
    out["young_mean"] = Y
    out["rs_mean"] = R
    eps = 1e-3
    out["log2FC_rs_vs_young"] = np.log2((out["rs_mean"]+eps)/(out["young_mean"]+eps))
    out.to_csv("/home/ibab/1_POJECTS/Telos_NGS_Project/1_data_raw/1_preprocessing_and_intermediate/atac_atlas_with_signal.csv",index=False)
    print("Wrote atac_atlas_with_signal.csv", out.shape)
#--------------------------------------------------------------------------------#

if __name__ == "__main__":
    main()
