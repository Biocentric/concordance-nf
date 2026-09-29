/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONCORDANCE - main analysis workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK               } from '../subworkflows/local/input_check'
include { BAM_CRAM_EXTRACT_SOMALIER } from '../subworkflows/local/bam_cram_extract_somalier'
include { SOMALIER_RELATE           } from '../modules/nf-core/somalier/relate/main'
include { SOMALIER_ANCESTRY         } from '../modules/nf-core/somalier/ancestry/main'
include { SOMALIER_CONTAMINATION    } from '../modules/local/somalier_contamination/main'
include { CONCORDANCE_REPORT        } from '../modules/local/concordance_report/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'

workflow CONCORDANCE {

    main:

    // ---- parameter validation -------------------------------------------
    if ( !params.input ) {
        error "--input is required (samplesheet CSV). Run with --help for the column spec."
    }
    if ( !params.fasta ) {
        error "--fasta is required. somalier extract always needs the reference, " +
              "and CRAM cannot be decoded without it."
    }

    def genome_cfg = params.genomes.containsKey(params.genome) ? params.genomes[params.genome] : null
    if ( !params.sites && !genome_cfg ) {
        error "Unknown --genome '${params.genome}' and no --sites given. " +
              "Known genomes: ${params.genomes.keySet().join(', ')}"
    }
    def sites_path = params.sites ?: genome_cfg.sites
    def fai_path   = params.fai   ?: "${params.fasta}.fai"

    if ( params.run_ancestry && !params.ancestry_labels ) {
        error "--run_ancestry needs --ancestry_labels (and --ancestry_labelled_dir)."
    }

    ch_samplesheet = Channel.value( file(params.input, checkIfExists: true) )
    ch_fasta       = Channel.value( [ [id: 'reference'],    file(params.fasta, checkIfExists: true) ] )
    ch_fai         = Channel.value( [ [id: 'reference'],    file(fai_path,     checkIfExists: true) ] )
    ch_sites       = Channel.value( [ [id: params.genome],  file(sites_path,   checkIfExists: true) ] )

    // ---- samplesheet -> sketches ----------------------------------------
    // a plain path, not the channel: INPUT_CHECK needs the parent directory
    INPUT_CHECK ( file(params.input, checkIfExists: true) )

    BAM_CRAM_EXTRACT_SOMALIER ( INPUT_CHECK.out.inputs, ch_fasta, ch_fai, ch_sites )

    // ---- historical sketch registry --------------------------------------
    // A sketch is ~200 KB and permanent. Comparing only within a batch is the
    // single biggest missed opportunity in routine identity QC: the swap you
    // care about is usually against a sample from three months ago.
    ch_registry = params.registry
        ? Channel.fromPath("${params.registry}/**.somalier", checkIfExists: true)
        : Channel.empty()

    ch_all_sketches = BAM_CRAM_EXTRACT_SOMALIER.out.sketches
        .map { _meta, sketch -> sketch }
        .mix( ch_registry )
        .collect()
        .map { sketches -> sketches.sort { a, b -> file(a).name <=> file(b).name } }

    // ---- all-vs-all relate ------------------------------------------------
    ch_ped    = params.ped           ? file(params.ped,           checkIfExists: true) : []
    ch_groups = params.sample_groups ? file(params.sample_groups, checkIfExists: true) : []

    SOMALIER_RELATE (
        ch_all_sketches.map { sketches -> [ [id: 'cohort'], sketches, ch_ped ] },
        ch_groups
    )

    // ---- contamination ----------------------------------------------------
    SOMALIER_CONTAMINATION (
        ch_all_sketches
            .filter { params.run_contamination }
            .map { sketches -> [ [id: 'cohort'], sketches ] },
        ch_sites
    )

    ch_contam_pairs   = SOMALIER_CONTAMINATION.out.pairs_tsv.map   { _m, f -> f }.ifEmpty( [] )
    ch_contam_samples = SOMALIER_CONTAMINATION.out.samples_tsv.map { _m, f -> f }.ifEmpty( [] )

    // ---- ancestry (optional) ----------------------------------------------
    ch_labelled = params.ancestry_labelled_dir
        ? Channel.fromPath("${params.ancestry_labelled_dir}/**.somalier", checkIfExists: true).collect()
        : Channel.value( [] )

    SOMALIER_ANCESTRY (
        ch_all_sketches
            .filter { params.run_ancestry }
            .map { sketches -> [ [id: 'cohort'], sketches ] },
        ch_labelled.map { labelled ->
            [ [id: 'labels'], params.ancestry_labels ? file(params.ancestry_labels, checkIfExists: true) : [], labelled ]
        }
    )

    ch_ancestry = SOMALIER_ANCESTRY.out.tsv.map { _m, f -> f }.ifEmpty( [] )

    // ---- verdicts against the declared identity structure ------------------
    ch_report_script = Channel.value( file("${projectDir}/bin/concordance_report.py", checkIfExists: true) )

    CONCORDANCE_REPORT (
        SOMALIER_RELATE.out.pairs_tsv.join( SOMALIER_RELATE.out.samples_tsv, failOnMismatch: true ),
        ch_samplesheet,
        ch_contam_pairs,
        ch_contam_samples,
        ch_ancestry,
        ch_report_script
    )

    // ---- MultiQC -----------------------------------------------------------
    ch_multiqc_files = Channel.empty()
        .mix( SOMALIER_RELATE.out.pairs_tsv.map   { _m, f -> f } )
        .mix( SOMALIER_RELATE.out.samples_tsv.map { _m, f -> f } )
        .mix( CONCORDANCE_REPORT.out.multiqc )
        .collect()
        .ifEmpty( [] )

    ch_multiqc_config = file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)

    MULTIQC (
        ch_multiqc_files
            .filter { !params.skip_multiqc }
            .map { files -> [ [id: 'concordance'], files, ch_multiqc_config, [], [], [] ] }
    )

    // ---- versions ----------------------------------------------------------
    channel.topic('versions')
        .map { _process, tool, version -> "  ${tool}: ${version}" }
        .unique()
        .collectFile(
            name:     'software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            sort:     true,
            seed:     "${workflow.manifest.name}:\n"
        )

    emit:
    pairs   = CONCORDANCE_REPORT.out.pairs
    samples = CONCORDANCE_REPORT.out.samples
    summary = CONCORDANCE_REPORT.out.summary
    html    = CONCORDANCE_REPORT.out.html
}
