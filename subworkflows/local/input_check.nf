//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_variant_channel(it) }
        .set { variants }

    emit:
    variants                                  // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_variant_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id             = row.sample

    // add path(s) of the fastq file(s) to the meta map
    def variants_meta = []
    if (!file(row.variants_file).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> variants_file file does not exist!\n${row.variants_file}"
    }
    if (!file(row.genome_ref_file).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> genome_ref_file file does not exist!\n${row.genome_ref_file}"
    }
    variants_meta = [ meta, file(row.variants_file), file(row.genome_ref_file) ]
    return variants_meta
}
