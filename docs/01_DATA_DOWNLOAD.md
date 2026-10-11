# Data and reference setup

Large files remain local under `1_data_raw/`; they are excluded from Git.

## Experimental source

- GEO accession recorded in the project notes: **GSE175533**.
- The local tree has processed ATAC bigWig files, atlas BED files, and spreadsheet/archives. These files are not included in this repository.
- Use the GEO record and the paper's sample metadata to verify each sample's cell condition and passage label. Do not rely only on filenames if the metadata differ.

## Reference genome

The current workflow uses UCSC hg38/GRCh38 chromosome names (`chr1`, `chr2`, …), consistent with the project BED files.

Expected local paths:

```text
1_data_raw/0_reference/hg38/hg38.fa
1_data_raw/0_reference/hg38/hg38.fa.fai
```

Download and index (requires `wget`, `gunzip`, and a working `samtools` installation):

```bash
mkdir -p 1_data_raw/0_reference/hg38
cd 1_data_raw/0_reference/hg38
wget -O hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -f hg38.fa.gz
samtools faidx hg38.fa
```

Verify:

```bash
ls -lh hg38.fa hg38.fa.fai
head -n 5 hg38.fa.fai
```

Use an already working Samtools environment if installing Samtools into the AlphaGenome environment causes dependency conflicts. Samtools only needs to index the FASTA here; it need not be installed into every Python environment.

## Directory creation

From the repository root, run:

```bash
bash scripts/create_local_dirs.sh
```

## Data governance

Document accession IDs, source URLs, retrieval date, genome build, sample mapping, and any processing steps. Do not commit raw GEO archives, bigWig tracks, genome FASTA/index files, or intermediates. Check the data provider's terms before redistributing any processed data.
