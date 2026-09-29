/*
 * CONCORDANCE_REPORT
 *
 * Turns somalier's raw relatedness numbers into verdicts against the identity
 * structure declared in the samplesheet. This is the piece neither somalier
 * nor nf-core provides: somalier tells you what the data says, not whether it
 * agrees with what you expected.
 *
 * The report script is staged as a regular input file rather than relying on
 * bin/ being on PATH. bin/ requires the execute bit, which does not survive on
 * every filesystem a pipeline gets checked out onto (notably FUSE mounts), and
 * a staged path works identically under every executor and container engine.
 */
process CONCORDANCE_REPORT {
    tag "${meta.id}"
    label 'process_single'

    // conda/container set in conf/modules.config

    input:
    tuple val(meta), path(pairs_tsv), path(samples_tsv)
    path samplesheet
    path contamination_pairs
    path contamination_samples
    path ancestry_tsv
    path report_script

    output:
    tuple val(meta), path("concordance.pairs.tsv"),    emit: pairs
    tuple val(meta), path("concordance.samples.tsv"),  emit: samples
    tuple val(meta), path("concordance.summary.json"), emit: summary
    tuple val(meta), path("concordance_report.html"),  emit: html
    path  "*_mqc.tsv",                                 emit: multiqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args              = task.ext.args ?: ''
    def contam_pairs_arg  = contamination_pairs   ? "--contamination-pairs ${contamination_pairs}"     : ''
    def contam_samp_arg   = contamination_samples ? "--contamination-samples ${contamination_samples}" : ''
    def ancestry_arg      = ancestry_tsv          ? "--ancestry ${ancestry_tsv}"                       : ''
    """
    python3 ${report_script} \\
        --pairs ${pairs_tsv} \\
        --samples ${samples_tsv} \\
        --samplesheet ${samplesheet} \\
        ${contam_pairs_arg} \\
        ${contam_samp_arg} \\
        ${ancestry_arg} \\
        --run-name '${meta.id}' \\
        --min-sites ${params.min_sites} \\
        --min-depth ${params.min_depth} \\
        --identity-min-relatedness ${params.identity_min_relatedness} \\
        --identity-max-ibs0-rate ${params.identity_max_ibs0_rate} \\
        --first-degree-min-relatedness ${params.first_degree_min_relatedness} \\
        --second-degree-min-relatedness ${params.second_degree_min_relatedness} \\
        --contamination-warn ${params.contamination_warn} \\
        --contamination-fail ${params.contamination_fail} \\
        ${params.fail_on_discordance ? '--fail-on-discordance' : ''} \\
        ${args} \\
        --outdir .
    """

    stub:
    // deliberately runs the real script: a stub run is the cheapest way to
    // prove the script stages and executes, and it copes with the empty TSVs
    // the upstream somalier stubs produce
    """
    python3 ${report_script} \\
        --pairs ${pairs_tsv} \\
        --samples ${samples_tsv} \\
        --samplesheet ${samplesheet} \\
        --run-name '${meta.id}' \\
        --outdir .
    """
}
