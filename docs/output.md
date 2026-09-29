# Output

```
results/
├── concordance/
│   ├── concordance_report.html      <- start here
│   ├── concordance.pairs.tsv
│   ├── concordance.samples.tsv
│   ├── concordance.summary.json
│   ├── concordance_pairs_mqc.tsv
│   └── concordance_samples_mqc.tsv
├── sketches/
│   └── <sample>.somalier            <- archive these
├── somalier/
│   ├── somalier.html                <- somalier's own interactive plots
│   ├── somalier.pairs.tsv
│   ├── somalier.samples.tsv
│   ├── somalier.contamination.samples.tsv
│   └── somalier.contamination.pairs.tsv
├── multiqc/
└── pipeline_info/
    ├── software_versions.yml
    ├── report.html, timeline.html, trace.txt
```

## concordance_report.html

Self-contained, no JavaScript, no external assets — safe to email or attach to a
ticket. Summary cards, then the pairs and samples needing attention, then a
relatedness matrix, then the full tables. The decision thresholds actually used
are printed in the footer.

## concordance.pairs.tsv

One row per evaluated pair. Every declared pair appears; registry pairs appear
only when they are not plainly unrelated, so a large registry does not drown the
table in N² uninteresting rows.

| Column | Notes |
|---|---|
| `sample_a`, `sample_b` | |
| `subject_a`, `subject_b` | from the samplesheet; blank for registry sketches |
| `expectation` | `same_subject`, `same_family`, `unrelated`, `unknown` |
| `relatedness` | somalier's KING-robust style estimate. ~1 identical, ~0.5 first-degree, ~0 unrelated |
| `ibs0`, `ibs0_rate` | hom-ref/hom-alt sites; impossible within one individual |
| `ibs2` | sites where both genotypes agree |
| `concordance` | somalier's own concordance. Reported, not used for the decision |
| `shared_hets`, `shared_hom_alts` | |
| `n_sites` | sites usable in *both* samples |
| `observed_class` | `IDENTICAL`, `PARENT_CHILD`, `SIBLING`, `SECOND_DEGREE`, `UNRELATED`, `INCONCLUSIVE` |
| `verdict` | observed class judged against the expectation |
| `severity` | `pass`, `warn`, `fail` |
| `note` | why, when it is not a pass |

## concordance.samples.tsv

One row per sample, including registry sketches (`in_samplesheet` = `no`).
Carries `n_sites`, `gt_depth_mean`, declared vs inferred sex with the evidence
used, `contamination_charr`, the top pairwise contaminator, ancestry when
enabled, and a `status` with `notes`.

A named `top_contaminator` is only reported when the sample's own CHARR is above
the warn threshold — below that, the anchored pairwise MLE is noise, and an
estimate pinned at 1.0 is the optimiser hitting its bound rather than a real
attribution. Both cases are suppressed or marked rather than presented as a
finding.

## concordance.summary.json

Machine-readable. Counts, the thresholds actually applied, `verdict_counts`,
`sample_status_counts`, the full `failures` and `sample_failures` lists, and
`missing_samples` — declared in the samplesheet but no sketch came back, which
usually means extraction failed or the sites file does not match the reference.

Gate a pipeline on it with `--fail_on_discordance`, which exits non-zero when
any verdict is `fail`.

## sketches/

The per-sample `.somalier` sketches, ~200 KB each. Archive them and pass the
directory as `--registry` on later runs. This is the single highest-value thing
in the output directory: the swap you care about is usually against a sample
from months ago, and a batch-scoped comparison can never find it.

They contain ~17k genotypes and are fully identifying. Treat them as personal
data.

## Interpreting a failure

**`SAMPLE_SWAP`** — two libraries declared as one subject disagree. Check the
per-sample rows first: if either has a low site count or low depth, the pair may
be underpowered rather than swapped. If both look healthy, it is a real
labelling error; `somalier.pairs.tsv` will usually show which *other* sample
each one actually matches.

**`UNEXPECTED_DUPLICATE`** — two declared-distinct subjects are genotypically
identical. Either the same material was sequenced twice under two identifiers,
or one library was mislabelled. Check `batch` and whether they are adjacent on a
plate.

**`UNEXPECTED_RELATEDNESS`** — relatedness between declared-unrelated subjects.
Before assuming an undocumented pedigree relationship, check contamination:
cross-contamination mimics relatedness, and the pairwise contamination table
distinguishes them because contamination is directional and relatedness is not.

**`INCONCLUSIVE`** — almost always depth or a contig-naming mismatch between the
sites file and the reference (`chr1` vs `1`), which yields zero usable sites
rather than an error.
