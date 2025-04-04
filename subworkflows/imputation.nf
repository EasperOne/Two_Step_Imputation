#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Include modules
include { CONVERT_REFP } from '../modules/convert_refp'
include { IMPUTE as IMPUTE_MINIMAC4 } from '../modules/impute'
include { IMPUTE_BEAGLE } from '../modules/impute_beagle'

workflow IMPUTATION {
    take:
    ch_phased_vcf_for_imputation // [meta, vcf, index] from phasing
    ch_input_vcfs                // Original input [meta, vcf] - for cohort information

    main:
    // Create a list of unique cohorts
    ch_cohorts = ch_input_vcfs
        .map { meta, vcf -> meta.cohort }
        .unique()
        .collect()
    
    if (params.imputer == "minimac4") {
        CONVERT_REFP(ch_phased_vcf_for_imputation)
        
        // Create a channel with reference panels by cohort and chromosome
        ch_ref_panels = CONVERT_REFP.out.reference
            .map { meta, panel -> [meta.cohort, meta.chr, meta, panel] }
        
        // Prepare phased VCFs for combination
        ch_phased_for_combo = ch_phased_vcf_for_imputation
            .map { meta, vcf, index -> [meta.cohort, meta.chr, meta, vcf, index] }
        
        // Create all valid combinations of phased VCFs and reference panels
        ch_impute_jobs = ch_phased_for_combo
            .combine(ch_ref_panels, by: 1)  // Join by chromosome first
            .filter { chr, phased_cohort, phased_meta, phased_vcf, phased_index, ref_cohort, ref_meta, ref_panel ->
                // Only impute when cohorts are different
                phased_cohort != ref_cohort
            }
            .map { chr, phased_cohort, phased_meta, phased_vcf, phased_index, ref_cohort, ref_meta, ref_panel ->
                // Create a new meta with target_cohort information required by post_imputation
                def newMeta = phased_meta.clone()
                newMeta.target_cohort = ref_cohort
                
                // Return format that matches IMPUTE_MINIMAC4 input: [meta, vcf, index, panel]
                [newMeta, phased_vcf, phased_index, ref_panel]
            }
        
        IMPUTE_MINIMAC4(ch_impute_jobs)
        ch_final_imputed_vcf = IMPUTE_MINIMAC4.out.imputed_vcf
    } 
    else if (params.imputer == "beagle") {
        // For Beagle imputation
        ch_phased_for_beagle = ch_phased_vcf_for_imputation
            .map { meta, vcf, index -> 
                // Clone and update meta to include required fields for POST_IMPUTATION
                def newMeta = meta.clone()
                newMeta.target_cohort = "${meta.cohort}_beagle"  // Using the same cohort as a pseudo-target since beagle doesn't use references
                return [newMeta, vcf, index]
            }
        
        IMPUTE_BEAGLE(ch_phased_for_beagle)
        
        // Ensure the output has the same structure as IMPUTE_MINIMAC4
        ch_final_imputed_vcf = IMPUTE_BEAGLE.out.imputed_vcf
            .combine(IMPUTE_BEAGLE.out.imputed_index, by: 0)
            .map { meta, vcf, index -> [meta, vcf, index] }
    }
    else {
        error "Unsupported imputer: ${params.imputer}. Choose 'minimac4' or 'beagle'."
    }

    emit:
    imputed_vcf = ch_final_imputed_vcf // [meta, vcf, index] format for POST_IMPUTATION
}
