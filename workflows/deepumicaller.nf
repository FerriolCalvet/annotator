/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowFgcons.initialise(params, log)

// TODO nf-core:
// Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.ref_fasta, params.targetsfile  ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.input) { ch_input = file(params.input) } else {
    exit 1, 'Input samplesheet not specified!'
    }

if (params.ref_fasta) {
    ch_ref_fasta = Channel.fromPath(params.ref_fasta).collect()

    // define additional fasta file names
    ch_ref_fasta_file = file(params.ref_fasta, checkIfExists: true)
    ch_ref_fasta_fai_index = file("${ch_ref_fasta_file}.fai", checkIfExists: true)
    ch_ref_fasta_dict = file("${ch_ref_fasta_file.parent/ch_ref_fasta_file.baseName}.dict", checkIfExists: true)

} else {
    log.error "No reference FASTA was specified (--ref_fasta)."
    exit 1
    }


// The index directory is the directory that contains the FASTA
ch_ref_index_dir = ch_ref_fasta.map { it -> it.parent }
// TODO
// check if the index file for the reference genome is present
// if (ch_ref_index_dir) { file("${file(params.ref_fasta).parent}/${file(params.ref_fasta).name}.amb", checkIfExists: true) }



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                                                      } from '../subworkflows/local/input_check'


// include { BBGPOSTANALYSIS                     as BBGPOSTANALYSIS               } from '../modules/local/BBGPOSTANALYSIS/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

// Download annotation cache if needed
include { PREPARE_CACHE                                                   } from '../subworkflows/local/prepare_cache/main'


// Annotation
include { VCF_ANNOTATE_ALL                  as VCFANNOTATELOW          } from '../subworkflows/local/vcf_annotate_all/main'
include { VCF_ANNOTATE_ALL                  as VCFANNOTATEMED          } from '../subworkflows/local/vcf_annotate_all/main'
include { VCF_ANNOTATE_ALL                  as VCFANNOTATEHIGH         } from '../subworkflows/local/vcf_annotate_all/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DEEPUMICALLER {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    // Download Ensembl VEP cache if needed
    // Assuming that if the cache is provided, the user has already downloaded it
    ensemblvep_info = params.vep_cache ? [] : Channel.of([ [ id:"${params.vep_genome}.${params.vep_cache_version}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
    if (params.download_cache) {
        PREPARE_CACHE(ensemblvep_info)
        vep_cache = PREPARE_CACHE.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }
        ch_versions = ch_versions.mix(PREPARE_CACHE.out.versions)
    } else {
        vep_cache = params.vep_cache
    }
    vep_extra_files = []
    
    
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


        VCFANNOTATEHIGH(CALLINGVARDICTHIGH.out.vcf,
                            ch_ref_fasta,
                            "GRCh38",
                            "homo_sapiens", 
                            "108",
                            vep_cache,
                            vep_extra_files)
        ch_versions = ch_versions.mix(VCFANNOTATEHIGH.out.versions.first())

        CALLINGVARDICTHIGH.out.vcf.map{it -> it[1]}.set { mutation_files_high }


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowFgcons.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/