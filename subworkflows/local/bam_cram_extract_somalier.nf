/*
 * BAM_CRAM_EXTRACT_SOMALIER
 *
 * The alignment-input counterpart to nf-core's `vcf_extract_relate_somalier`,
 * which handles VCFs only and requires a PED. Modelled on that subworkflow's
 * index-branching, extended to cover what the identity use case actually needs:
 *
 *   BAM/CRAM with an index   -> straight to somalier extract
 *   BAM/CRAM without one     -> SAMTOOLS_INDEX first
 *   VCF already bgzipped+idx -> straight through
 *   VCF otherwise            -> HTSLIB_BGZIPTABIX (compress and/or index)
 *   BCF                      -> index required up front (see input_check)
 *
 * Extracting from the alignment is the preferred path. Sketching a jointly
 * called VCF inherits the information sharing that joint calling introduces
 * across samples, which inflates apparent concordance.
 */

include { SAMTOOLS_INDEX    } from '../../modules/nf-core/samtools/index/main'
include { HTSLIB_BGZIPTABIX } from '../../modules/nf-core/htslib/bgziptabix/main'
include { SOMALIER_EXTRACT  } from '../../modules/nf-core/somalier/extract/main'

workflow BAM_CRAM_EXTRACT_SOMALIER {
    take:
    ch_inputs    // channel: [ val(meta), path(file), path(index) | [] ]
    ch_fasta     // channel: [ val(meta2), path(fasta) ]
    ch_fai       // channel: [ val(meta3), path(fai) ]
    ch_sites     // channel: [ val(meta4), path(sites_vcf) ]

    main:

    ch_branch = ch_inputs
        .branch { meta, _datafile, index ->
            // order matters: first match wins, last arm must be unconditional
            aln_indexed:  (meta.filetype == 'bam' || meta.filetype == 'cram') && index
            aln_bare:     (meta.filetype == 'bam' || meta.filetype == 'cram')
            var_indexed:  index as boolean
            var_bare:     true
        }

    // ---- alignments needing an index ------------------------------------
    ch_aln_bare = ch_branch.aln_bare.map { meta, datafile, _index -> [ meta, datafile ] }

    SAMTOOLS_INDEX ( ch_aln_bare )

    ch_aln_freshly_indexed = ch_aln_bare
        .join( SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true )

    // ---- variants needing compression and/or an index --------------------
    HTSLIB_BGZIPTABIX (
        ch_branch.var_bare.map { meta, datafile, _index -> [ meta, datafile, [], [] ] },
        'compress',
        true,
        'vcf'
    )

    ch_var_prepared = HTSLIB_BGZIPTABIX.out.output
        .join( HTSLIB_BGZIPTABIX.out.index, failOnMismatch: true, failOnDuplicate: true )

    // ---- one stream into extract ----------------------------------------
    ch_for_extract = ch_branch.aln_indexed
        .mix( ch_aln_freshly_indexed )
        .mix( ch_branch.var_indexed )
        .mix( ch_var_prepared )

    SOMALIER_EXTRACT ( ch_for_extract, ch_fasta, ch_fai, ch_sites )

    emit:
    sketches = SOMALIER_EXTRACT.out.extract   // [ meta, *.somalier ]
    staged   = ch_for_extract                 // [ meta, file, index ]
}
