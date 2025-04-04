#!/usr/bin/env nextflow

// Include modules
include { ALIGNMENT } from '../modules/alignment'
include { VCF_TO_BCF } from '../modules/vcf_to_bcf'

workflow PRE_PHASE {
    take:
    ch_input_vcfs // Channel of [meta, vcf_file]

    main:
    // Alignment step
    ALIGNMENT(ch_input_vcfs)
    
    // Convert VCF to BCF
    VCF_TO_BCF(ALIGNMENT.out.vcf)

    emit:
    bcf = VCF_TO_BCF.out.bcf // Output [meta, bcf_file, bcf_index]
}
