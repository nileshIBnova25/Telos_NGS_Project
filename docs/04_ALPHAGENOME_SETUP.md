# AlphaGenome setup and safe API use

## Python environment

A starter Conda specification is provided in `envs/alphagenome-env.yml`. It is a convenience starting point, not a fully locked reproducibility file. After the environment is working, export and review a platform-appropriate lock/export file and document the tested package versions.

```bash
conda env create -f envs/alphagenome-env.yml
conda activate alphagenome-env
python -c "import alphagenome, pyfaidx; print('imports OK')"
```

If the environment already exists, update deliberately rather than blindly applying the file. For genome indexing, use a separate known-working Samtools environment if that avoids dependency conflicts.

## Store the key in the shell environment

Use a key obtained through the official AlphaGenome access flow and check current usage terms before making large numbers of requests. Never paste the key into source code, commit it, place it in a notebook, or print it to logs.

```bash
read -rsp "AlphaGenome API key: " ALPHA_GENOME_API_KEY
echo
export ALPHA_GENOME_API_KEY
```

This variable lasts only for the current shell session. The script may also accept `ALPHAGENOME_API_KEY` if that is the name used by the current implementation. Keep variable names consistent and document which one the script reads.

A `.env.example` file can document the variable name only; it must contain no real key. Local `.env` files are ignored.

## Track selection

- Query official AlphaGenome output metadata for human ATAC tracks.
- Search and inspect track metadata rather than assuming WI-38 exists.
- Prefer the best-matched biosample only if metadata substantiate the match; otherwise clearly label the proxy.
- Save track metadata alongside prediction scores so the analysis can be audited.
- A successful API response confirms a request worked; it does not establish biological validity.

## Prediction outputs

Save sequences, predictions, metadata, plots, logs and completion markers locally under `1_data_raw/3_alphagenome/`. Do not commit them. For each run, document the input interval (hg38, coordinate convention), sequence length, requested output, selected ontology/track, package version, date, and failure/retry behavior.

## Protect credentials

Before committing, review all files for accidental secrets. If a key was pasted into a chat, issue tracker, terminal transcript, or file shared with others, treat it as potentially exposed and revoke/rotate it if it is real. Removing it from the current file does not remove it from Git history if it has already been committed.
