/*
 * Parse and validate the samplesheet, and classify each row by input type.
 *
 * Uses Nextflow's native splitCsv rather than the nf-schema plugin, matching
 * the house style in gvanno-nf: it keeps the dependency surface small enough
 * that the pipeline runs offline without a plugin download.
 * assets/samplesheet_schema.json documents the same contract for editors and
 * for anyone who does want to validate externally.
 *
 * Note for Nextflow >= 25: the strict parser allows only process, workflow and
 * function declarations at script level, so constants live inside functions.
 */

def classify(String path) {
    def alignment_ext = ['.bam', '.cram', '.sam']
    def variant_ext   = ['.vcf', '.vcf.gz', '.gvcf', '.gvcf.gz', '.g.vcf', '.g.vcf.gz']
    def lower = path.toLowerCase()

    if ( alignment_ext.any { lower.endsWith(it) } ) {
        return lower.endsWith('.cram') ? 'cram' : 'bam'
    }
    if ( lower.endsWith('.bcf') || lower.endsWith('.bcf.gz') ) {
        return 'bcf'
    }
    if ( variant_ext.any { lower.endsWith(it) } ) {
        return 'vcf'
    }
    return null
}

def normalise_sex(String value) {
    def v = (value ?: '').trim().toLowerCase()
    if ( v in ['m', 'male', '1'] )                  return 'male'
    if ( v in ['f', 'female', '2'] )                return 'female'
    if ( v in ['', 'u', 'unknown', '0', '-9', 'na'] ) return 'unknown'
    error "Samplesheet: unrecognised sex '${value}'. Use male, female or unknown."
}

/*
 * Resolve a path from the samplesheet. Absolute paths and URLs are taken as
 * given; a relative path is tried against the launch directory first (the
 * Nextflow default) and then against the samplesheet's own directory, so a
 * samplesheet that sits next to its data stays portable when run from elsewhere.
 */
def resolve_path(String raw, sheet_dir) {
    def candidate = file(raw)
    if ( candidate.exists() ) return candidate
    def beside = sheet_dir.resolve(raw)
    if ( beside.exists() ) return beside
    error "File not found: '${raw}' (tried the launch directory and ${sheet_dir})"
}

workflow INPUT_CHECK {
    take:
    samplesheet   // a plain path, NOT a channel - we need its parent directory
                  // eagerly to resolve samplesheet-relative entries

    main:

    def sheet_file = file(samplesheet, checkIfExists: true)
    def sheet_dir  = sheet_file.parent

    ch_rows = Channel
        .fromPath(sheet_file, checkIfExists: true)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->

            if ( !row.sample )    error "Samplesheet row is missing a 'sample' value: ${row}"
            if ( !row.subject )   error "Samplesheet row '${row.sample}' is missing a 'subject' value. " +
                                        "subject is what the concordance check is tested against - " +
                                        "use the sample id itself if every row is a distinct individual."
            if ( !row.alignment ) error "Samplesheet row '${row.sample}' is missing an 'alignment' value."

            def input_file = resolve_path(row.alignment as String, sheet_dir)
            def filetype   = classify(input_file.name)
            if ( !filetype ) {
                error "Samplesheet row '${row.sample}': unsupported input '${input_file.name}'. " +
                      "Expected BAM, CRAM, SAM, VCF, BCF or GVCF."
            }

            def index_file = row.index ? resolve_path(row.index as String, sheet_dir) : []

            // BCF is already binary and compressed, so there is nothing for
            // bgzip to do and no safe way to derive an index from the
            // extension alone. Ask for one rather than guessing.
            if ( filetype == 'bcf' && !index_file ) {
                error "Samplesheet row '${row.sample}': BCF input needs an 'index' column value " +
                      "(bcftools index ${input_file.name})."
            }

            def meta = [
                id:       row.sample,
                subject:  row.subject,
                family:   row.family ?: row.subject,
                sex:      normalise_sex(row.sex as String),
                batch:    row.batch  ?: 'default',
                filetype: filetype
            ]
            tuple(meta, input_file, index_file)
        }

    // Fail fast on duplicate sample ids: somalier keys sketches by the sample
    // name embedded in the file, so a duplicate silently overwrites rather than
    // erroring, and the resulting matrix is then quietly wrong.
    ch_rows
        .map { meta, _f, _i -> meta.id }
        .toList()
        .map { ids ->
            def dupes = ids.countBy { it }.findAll { _k, v -> v > 1 }.keySet()
            if ( dupes ) {
                error "Samplesheet contains duplicate sample ids: ${dupes.join(', ')}"
            }
            ids.size()
        }
        .subscribe { n -> log.info "[concordance-nf] samplesheet: ${n} samples" }

    emit:
    inputs = ch_rows     // [ meta, file, index|[] ]
}
