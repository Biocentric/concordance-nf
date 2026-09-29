# Stub-test placeholders

The `.cram`, `.crai` and `.bam` files here are **deliberately empty**. They exist
so `-profile test_stub` can exercise the BAM/CRAM input branches — indexed,
unindexed, and mixed with a VCF — under `-stub`, where no tool ever opens them.

They are not valid alignments and will fail any real run.

For a genuine end-to-end test against real somalier, use `-profile test`, which
reads the synthetic cohort in [`../test_data/`](../test_data/).
