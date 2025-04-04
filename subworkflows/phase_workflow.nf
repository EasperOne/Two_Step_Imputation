#!/usr/bin/env nextflow

// Include modules
include { PHASE as PHASE_EAGLE } from '../modules/phase'
include { PHASE_SHAPEIT5 } from '../modules/phase_shapeit5'
include { BCF_TO_VCF } from '../modules/bcf_to_vcf'

workflow PHASE_WORKFLOW {
    take:
    ch_bcf // Channel of [meta, bcf_file, bcf_index]

    main:
    ch_phased_vcf = Channel.empty()

    // --- Phasing Step (Conditional) ---
    if (params.phaser == 'eagle') {
        PHASE_EAGLE(ch_bcf)
        ch_phased_vcf = PHASE_EAGLE.out.phased_vcf // Eagle outputs VCF.gz
    } else if (params.phaser == 'shapeit5') {
        PHASE_SHAPEIT5(ch_bcf)
        BCF_TO_VCF(PHASE_SHAPEIT5.out.phased_bcf) // Convert SHAPEIT5 BCF to VCF
        ch_phased_vcf = BCF_TO_VCF.out.converted_vcf
    } else {
        error "Invalid phaser specified: ${params.phaser}. Choose 'eagle' or 'shapeit5'."
    }

    ch_phased_vcf.view { "Phased VCF ready for Imputation/RefPanel: ${it[0].id}" }

    emit:
    phased_vcf = ch_phased_vcf // Output [meta, vcf, index]
}
