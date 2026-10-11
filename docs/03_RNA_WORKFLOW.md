# RNA-seq workflow — planned, not yet documented as completed

RNA-seq integration is a later phase. The project notes describe the intended sequence:

1. Verify the GSE175533 RNA sample metadata and map samples to the same biological conditions where appropriate.
2. Perform RNA-seq QC and preprocessing, followed by alignment or quantification and gene-level counts/normalization appropriate to the design.
3. Run differential-expression analysis with a documented model and contrast. Record biological replicate counts and potential confounders.
4. Request AlphaGenome RNA-related outputs only after identifying the relevant output type, biosample/track context, sequence window, and aggregation method.
5. Compare prediction to experimental measurements using predeclared metrics and analyze results by region/gene class.
6. If linking ATAC peaks to genes, document the association method and distance/window rules. Accessibility changes and RNA expression changes are associations, not automatically proof of enhancer–gene causality.
7. Add motif, transcription-factor, or pathway analysis only after the underlying differential and association results have passed QC.

Do not describe this workflow as complete until the inputs, scripts, outputs and checks are available and documented.
