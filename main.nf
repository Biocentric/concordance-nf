#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Biocentric/concordance-nf
    https://github.com/Biocentric/concordance-nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SNP concordance / sample-identity QC for WGS, WES and targeted data.

    Wraps the nf-core somalier modules with the BAM/CRAM wiring they lack
    (nf-core ships only a VCF-input subworkflow), a samplesheet that carries
    the *expected* identity structure, and a report that turns somalier's
    raw relatedness numbers into per-pair verdicts against that expectation.

    somalier is the work of Brent Pedersen and Aaron Quinlan.
    https://github.com/brentp/somalier - doi:10.1186/s13073-020-00761-2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { CONCORDANCE } from './workflows/concordance'

def helpMessage() {
    log.info """
    concordance-nf - SNP concordance and sample-identity QC (somalier)

    Usage:
        nextflow run Biocentric/concordance-nf -profile docker \\
            --input samplesheet.csv \\
            --fasta /ref/GRCh38.fa \\
            --genome GRCh38 \\
            --outdir results

    Required:
        --input             Samplesheet CSV. Columns:
                              sample     (required) unique library/run id
                              subject    (required) the individual this sample
                                         should belong to - rows sharing a
                                         subject MUST come back identical
                              alignment  (required) BAM / CRAM / VCF / BCF / GVCF
                              index      (optional) .bai/.crai/.csi/.tbi;
                                         generated on the fly when omitted
                              family     (optional) relatedness between
                                         different subjects in one family is
                                         expected, not flagged
                              sex        (optional) male | female | unknown
                              batch      (optional) free-text grouping label
        --fasta             Reference FASTA. Required for CRAM.

    Reference:
        --genome            GRCh38 | GRCh38_nochr | GRCh37 | hg19 | CHM13 | GRCh38_rna
                            (default: GRCh38 - UCSC-style `chr` contigs)
        --fai               Defaults to <fasta>.fai
        --sites             Override the per-genome somalier sites VCF

    Sketch registry:
        --registry          Directory of previously extracted *.somalier files.
                            Mixed into the all-vs-all so new samples are checked
                            against project history, not just their own batch.
        --publish_sketches  Publish this run's sketches (default: true)

    Analysis:
        --run_contamination     somalier contamination, pairwise (default: true)
        --run_ancestry          somalier ancestry (default: false)
        --ancestry_labels       Labels TSV, required with --run_ancestry
        --ancestry_labelled_dir Directory of labelled *.somalier sketches
        --ped                   PED file for pedigree-aware relate
        --sample_groups         somalier -g groups file
        --skip_multiqc          Skip the MultiQC aggregation

    Decision thresholds:
        --min_sites                      Min usable sites per sample (1000)
        --min_depth                      Min mean genotype depth (5.0)
        --identity_min_relatedness       Identical if relatedness >= this (0.85)
        --identity_max_ibs0_rate         ...and IBS0 rate <= this (0.005)
        --first_degree_min_relatedness   (0.35)
        --second_degree_min_relatedness  (0.15)
        --contamination_warn             CHARR warn threshold (0.02)
        --contamination_fail             CHARR fail threshold (0.05)
        --fail_on_discordance            Exit non-zero on any FAIL (false)

    Tool version:
        --somalier_version    (default 0.3.5 - the nf-core modules pin 0.2.19,
                               which predates CHARR contamination and two
                               correctness fixes; see docs/usage.md)

    Other:
        --outdir            Default ./results
    """
}

workflow {
    if ( params.help ) {
        helpMessage()
        return
    }
    if ( params.version ) {
        log.info "${workflow.manifest.name} ${workflow.manifest.version}"
        return
    }
    CONCORDANCE()
}
