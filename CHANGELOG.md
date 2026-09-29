# Changelog

## 0.1.0 - 2026-09-29

Initial release.

- somalier-based SNP concordance / sample-identity QC for BAM, CRAM, VCF, BCF
  and GVCF input.
- Samplesheet carries the expected identity structure (`subject`, `family`), so
  the report produces verdicts rather than raw relatedness numbers.
- `BAM_CRAM_EXTRACT_SOMALIER`: the alignment-input subworkflow nf-core does not
  ship, including on-the-fly indexing.
- Persistent sketch registry via `--registry`, so new samples are compared
  against project history rather than only their own batch.
- Pinned to somalier 0.3.5 rather than the 0.2.19 the nf-core modules carry.
  0.3.x brings CHARR contamination, the pairwise `contamination` subcommand, a
  long-read allele-counting fix (0.3.4) and a hom-ref/hom-alt counting fix
  (0.3.5), and renames `homozygous_concordance` to `concordance`.
- `SOMALIER_CONTAMINATION` local module: no nf-core module exists for the
  subcommand.
- Self-contained HTML report, per-pair and per-sample TSVs, machine-readable
  JSON summary, and MultiQC custom-content sections.
