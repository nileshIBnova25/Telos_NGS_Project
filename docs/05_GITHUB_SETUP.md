# Fresh Git and GitHub setup for TELOS

**Starting point:** `/home/ibab/1_POJECTS/Telos_NGS_Project` is a normal project directory and is **not yet a Git repository**. These steps initialize Git from scratch. Do not add a remote or push until the first commit has been reviewed.

## 1. Review files and secrets before tracking anything

```bash
cd /home/ibab/1_POJECTS/Telos_NGS_Project
find . -maxdepth 3 -type f | sort
```

Do not commit `1_data_raw/`, `Research_paper/`, local ODT notes, or `2_scripts/other_script_utils_used_all.txt`. Review every script for hard-coded API keys, tokens, passwords, private paths, and machine-specific configuration. If a real API key has been exposed, revoke/rotate it before publishing.

## 2. Copy the contributor starter files into the project

Extract the starter ZIP to a temporary location, then from the project root copy its files. Merge existing documents rather than silently overwriting them.

```bash
# Set this to the extracted starter folder on your machine.
STARTER=/path/to/telos_github_starter

cp -n "$STARTER/.gitignore" .
cp -n "$STARTER/README.md" .
cp -n "$STARTER/CONTRIBUTING.md" .
cp -n "$STARTER/.env.example" .
mkdir -p docs envs scripts
cp -rn "$STARTER/docs/." docs/
cp -rn "$STARTER/envs/." envs/
cp -n "$STARTER/scripts/create_local_dirs.sh" scripts/
chmod +x scripts/create_local_dirs.sh
```

If you used the suggested local download path, the ZIP may be in `~/Downloads/`. Adjust the `STARTER` path to wherever you extracted it.

## 3. Initialize a brand-new local Git repository

Run from the project root:

```bash
git init -b main
```

This creates the local Git repository. No GitHub remote exists yet, and nothing has been uploaded.

Check the current state:

```bash
git status --short --branch
```

## 4. Verify raw data and secrets are ignored

```bash
git check-ignore -v \
  1_data_raw/0_reference/hg38/hg38.fa \
  1_data_raw/0_reference/hg38/hg38.fa.fai \
  1_data_raw/1_preprocessing_and_intermediate/atac_atlas_with_signal.csv \
  1_data_raw/GSE175533_RAW/GSM5341071_PDL20_TP1_omni_A_S4_R1_001.filter.rmdup.bw \
  2_scripts/other_script_utils_used_all.txt
```

Each path should display a matching `.gitignore` rule. If it does not, stop and inspect `.gitignore` before staging. `.gitignore` does not untrack files that have already been staged or committed, so avoid staging data in the first place.

## 5. Stage only documentation and reviewed source code

Do not run `git add .` for the initial commit. Start with known, reviewed paths:

```bash
git add .gitignore README.md CONTRIBUTING.md .env.example
git add docs envs scripts
# Inspect these scripts for secrets before staging them:
git add 2_scripts/1_Preprocessing

git status --short
git diff --cached --stat
git diff --cached --check
git diff --cached
```

Check the full staged diff. It must not contain the FASTA, bigWigs, GEO tar archives, generated benchmark files, CSV/XLSX outputs, logs, downloaded papers, ODT notes, API keys, or private scratchpad. If anything unexpected appears, unstage it and fix the ignore rules before committing.

## 6. Configure Git identity if needed, then commit

Only configure this if Git reports that `user.name` or `user.email` is missing:

```bash
git config user.name "Your Name"
git config user.email "you@example.com"
```

Use the identity you want attached to project commits. Then make the first commit:

```bash
git commit -m "Initial commit: document TELOS analysis workflow"
```

## 7. Create the GitHub repository and push

On GitHub, create a **new empty repository** named `Telos_NGS_Project`. Keep it private while auditing content. Because the local repository already has its first commit, do not ask GitHub to initialize the remote with another README, `.gitignore`, or license.

Replace `YOUR-USERNAME` with the actual GitHub username:

```bash
git remote add origin https://github.com/YOUR-USERNAME/Telos_NGS_Project.git
git push -u origin main
```

If GitHub asks you to authenticate, use the supported browser or credential-manager flow; do not put a password or personal access token directly in a command or source file.

## 8. Verify the upload

On GitHub, confirm that the repository contains source scripts, guides, and environment configuration only. Confirm that `1_data_raw/`, `Research_paper/`, working ODT files, logs, API keys, and `other_script_utils_used_all.txt` are absent. Keep the repository private until you finish this audit and choose an appropriate code license.

### If a secret or large file was committed accidentally

Do not push. Revoke exposed credentials immediately. Removing a file in a later commit does not remove it from earlier Git history; the history must be cleaned before publishing.
