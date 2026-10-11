# Project aim and scientific scope

## Research question

How accurately can AlphaGenome predict experimentally observed ATAC-seq and later RNA-seq regulatory activity in WI-38 fibroblasts, and does accuracy vary in regions associated with replicative senescence?

## Hypothesis boundary

The comparison is intended to test which regulatory patterns are predictable from DNA sequence and which observed differences depend on cellular/epigenetic state. Young and late-passage WI-38 samples generally share the same reference sequence; passage-associated chromatin changes are not encoded as a separate DNA sequence provided to the model.

Do not describe this as direct prediction of the senescent condition. State what the model receives, which output track was requested, which experimental group is the ground truth, and which metric was used.

## Conditions currently recorded

- Young/proliferating: PDL20, three replicates (A/B/C).
- Late/replicative-senescence group: PDL46 plus PDL50, four replicates (two each).

These are the project's current operational group definitions. Any change must be documented and propagated consistently through the benchmark and its labels.

## Intended outputs

1. A reproducible set of benchmark genomic intervals.
2. Experimental young and late/RS ATAC measurements for these intervals.
3. AlphaGenome predictions from documented output tracks and biosample context.
4. Per-region and stratified comparisons, including failure-case analysis.
5. Later, RNA/TSS analyses only after the ATAC evaluation is validated.

## Interpretation risks to report

- A track from a generic fibroblast or another cell type is not equivalent to a WI-38, passage-matched experiment.
- Shuffled regions outside the atlas may be easier controls than biological gained/lost/stable peaks.
- A sequence-only prediction should not be interpreted as having observed passage number or senescence state.
- Reference build, chromosome naming, BED coordinate conventions, track units, normalization, and region-selection procedures must be kept consistent.
