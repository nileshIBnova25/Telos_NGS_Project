#!/usr/bin/env python3
"""Fetch AlphaGenome output metadata and identify fibroblast-related tracks."""
from pathlib import Path
import os
import re
import sys
import pandas as pd
from alphagenome.models import dna_client

ROOT = Path.cwd()
OUT_DIR = ROOT / "1_data_raw/4_rnaseq/outputs/alphagenome"


def main():
    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        raise RuntimeError(
            "Set ALPHAGENOME_API_KEY in this terminal first. Do not put the key in the script or commit it to Git."
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    model = dna_client.create(key)
    meta_obj = model.output_metadata(dna_client.Organism.HOMO_SAPIENS)
    meta = meta_obj.concatenate()
    all_path = OUT_DIR / "alphagenome_human_output_metadata.csv"
    meta.to_csv(all_path, index=False)
    print(f"Saved full model output metadata: {all_path}")
    print("Metadata columns:", list(meta.columns))

    low = meta.astype(str).apply(lambda col: col.str.lower())
    row_text = low.agg(" ".join, axis=1)
    # First search for fibroblast; if there are none, return related tissue/cell terms
    # rather than silently substituting an unrelated context.
    fibro = meta[row_text.str.contains("fibroblast", regex=False)].copy()
    if not fibro.empty:
        selected = fibro
        print(f"\nFound {len(fibro)} metadata rows containing 'fibroblast':")
    else:
        related = row_text.str.contains(r"mesenchymal|connective tissue|lung|skin", regex=True)
        selected = meta[related].copy()
        print("\nNo track row literally contains 'fibroblast'. Showing related metadata matches instead.")
        print("Do NOT treat these as WI-38-specific; review biosample names and ontology IDs.")

    related_path = OUT_DIR / "alphagenome_fibroblast_related_tracks.csv"
    selected.to_csv(related_path, index=False)
    print(f"Saved relevant metadata subset: {related_path}")
    if selected.empty:
        print("No fibroblast-related or fallback metadata rows found. Inspect the full metadata CSV.")
    else:
        # Show the biologically useful fields if available, otherwise print all columns.
        wanted = [c for c in ["output_type", "name", "ontology_curie", "biosample_name", "biosample_type", "strand", "data_source"] if c in selected.columns]
        print(selected[wanted or list(selected.columns)].drop_duplicates().head(100).to_string(index=False))
        terms = [c for c in ["ontology_curie", "biosample_name"] if c in selected.columns]
        if "ontology_curie" in selected.columns:
            term_table = selected[terms].drop_duplicates().sort_values(terms)
            term_path = OUT_DIR / "alphagenome_candidate_ontology_terms.csv"
            term_table.to_csv(term_path, index=False)
            print(f"\nOntology term shortlist written to: {term_path}")
    print("\nNext: choose only ontology CURIEs whose biosample descriptions are biologically appropriate. The model may not have WI-38 or senescent-fibroblast tracks.")

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
