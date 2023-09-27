//
// ANNOTATION
//

include { VCF_ANNOTATE_ENSEMBLVEP                       } from '../../nf-core/vcf_annotate_ensemblvep/main.nf'
// include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_MERGE } from '../../nf-core/vcf_annotate_ensemblvep/main.nf'
// include { VCF_ANNOTATE_SNPEFF                           } from '../../nf-core/vcf_annotate_snpeff/main.nf'

workflow VCF_ANNOTATE_ALL {
    take:
    tab          // channel: [ val(meta), tab, fasta ]
    vep_genome
    vep_species
    vep_cache_version
    vep_cache
    vep_extra_files

    main:
    reports = Channel.empty()
    vcf_ann = Channel.empty()
    tab_ann = Channel.empty()
    json_ann = Channel.empty()
    versions = Channel.empty()

    VCF_ANNOTATE_ENSEMBLVEP(tab, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

    reports = reports.mix(VCF_ANNOTATE_ENSEMBLVEP.out.reports)
    // vcf_ann = vcf_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi)
    vcf_ann = vcf_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf)
    tab_ann = tab_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.tab)
    json_ann = json_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.json)
    versions = versions.mix(VCF_ANNOTATE_ENSEMBLVEP.out.versions)
    // }

    emit:
    vcf_ann      // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab_ann
    json_ann
    reports      //    path: *.html
    versions     //    path: versions.yml
}