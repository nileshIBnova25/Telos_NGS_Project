# Contributing to TELOS

Thank you for helping make this analysis reproducible. The repository should contain code, documentation, environment specifications, and small non-sensitive tests—not the research data or generated outputs.

## Before making changes

1. Read `README.md` and the relevant workflow document under `docs/`.
2. Check the current branch and working tree:

   ```bash
   git status --short --branch
   git branch --show-current
   ```

3. Create a focused branch for a change:

   ```bash
   git switch -c docs/describe-step4a
   ```

## Data and secrets rules

- Do not add files from `1_data_raw/`.
- Do not commit API keys, `.env` files, credentials, access tokens, or screenshots/logs that might expose them.
- `2_scripts/other_script_utils_used_all.txt` is intentionally excluded: it is a local scratchpad and may contain credentials or unreviewed notes.
- Keep downloaded papers and private working documents local unless redistribution rights have been checked.
- Use public accession IDs and documented download steps instead of copying large files into the repository.

## Code expectations

- Keep scripts runnable from a documented working directory or derive paths from the script/repository location rather than relying on a single contributor's absolute `/home/...` path.
- Put reusable logging or path utilities in a shared module when appropriate.
- Make scripts fail clearly when inputs or required tools are missing.
- Record important parameters, package versions, random seeds, genome build, data accessions, and output schema.
- Write generated files under `1_data_raw/` so they remain outside version control.
- Do not mark a workflow step complete unless the output exists and has been checked.
- For API workflows, support environment-variable credentials and never print the key.

## Validation before committing

Run checks relevant to the change, then inspect the staged files:

```bash
git diff --check
git status --short
git diff --cached --stat
git diff --cached
```

For a first commit, also verify that large data and sensitive notes are ignored:

```bash
git check-ignore -v \
  1_data_raw/0_reference/hg38/hg38.fa \
  1_data_raw/GSE175533_RAW/GSM5341071_PDL20_TP1_omni_A_S4_R1_001.filter.rmdup.bw \
  2_scripts/other_script_utils_used_all.txt
```

Use clear commit messages, for example:

```text
docs: document ATAC benchmark construction
fix: validate benchmark BED coordinates
feat: add resumable AlphaGenome pilot
```

## Pull requests

Describe the question addressed, files changed, commands run, outputs checked, and known limitations. Do not attach raw data or keys. For scientific changes, explain how the modification affects the benchmark or interpretation.
