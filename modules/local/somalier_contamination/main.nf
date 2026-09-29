/*
 * somalier contamination  (requires somalier >= 0.3.2)
 *
 * There is no nf-core module for this subcommand - the nf-core somalier
 * modules pin 0.2.19, which predates it entirely.
 *
 * Two outputs:
 *   PREFIX.samples.tsv  per-sample CHARR estimate
 *                       #sample_name  n_sites_usable  contamination_charr
 *   PREFIX.pairs.tsv    Conpair-style anchored pairwise MLE, i.e. *which*
 *                       sample is contaminating which
 *                       #sample_name  anchor_sample  n_sites_usable  contamination_mle
 *
 * The pairwise estimate is what disambiguates a genuine unexpected-relatedness
 * finding from cross-contamination between two libraries: contamination shows
 * up as directional, relatedness does not.
 *
 * --sites must be the same site set used for extract; somalier validates the
 * site count in each sketch header and aborts on a mismatch.
 */
process SOMALIER_CONTAMINATION {
    tag "${meta.id}"
    label 'process_low'

    // conda/container are set centrally in conf/modules.config so the somalier
    // version is pinned in exactly one place across all four somalier processes

    input:
    tuple val(meta),  path(extracted, stageAs: "extracted/*")
    tuple val(meta2), path(sites)

    output:
    tuple val(meta), path("*.samples.tsv"), emit: samples_tsv
    tuple val(meta), path("*.pairs.tsv"),   emit: pairs_tsv
    tuple val("${task.process}"), val('somalier'), eval('somalier 2>&1 | sed -n \'s/.*version: \\([0-9.]*\\).*/\\1/p\''), emit: versions_somalier, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    somalier contamination \\
        --sites ${sites} \\
        -o ${prefix} \\
        ${args} \\
        extracted/*.somalier
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '#sample_name\\tn_sites_usable\\tcontamination_charr\\n' > ${prefix}.samples.tsv
    printf '#sample_name\\tanchor_sample\\tn_sites_usable\\tcontamination_mle\\n' > ${prefix}.pairs.tsv
    """
}
