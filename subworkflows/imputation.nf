#!/usr/bin/env nextflow

// Include modules
include { CONVERT_REFP } from '../modules/convert_refp'
include { IMPUTE as IMPUTE_MINIMAC4 } from '../modules/impute'
include { IMPUTE_BEAGLE } from '../modules/impute_beagle'

workflow IMPUTATION {
    take:
    ch_phased_vcf_for_imputation // [meta, vcf, index] from phasing
    ch_input_vcfs                // Original input [meta, vcf] - for cohort information

    main:
    // --- Imputation Setup (Part 1: Reference Panels) ---
    // Get cohort names from input VCFs
    ch_cohort_names = ch_input_vcfs
        .map { meta, _vcf -> meta.cohort }
        .unique()
        .collect()
        .view { "Cohort Names: $it" }

    // Create imputation pairs for cross-imputation
    ch_imputation_pairs = ch_cohort_names
        .flatMap { cohorts ->
            cohorts.combinations().findAll { it.size() == 2 }.collectMany { 
                def a = it[0]
                def b = it[1]
                // Create pairs in both directions [A,B] and [B,A]
                return [[a, b], [b, a]]
            }
        }
        .view { "Imputation Pair: Source ${it[0]}, Target ${it[1]}" }

    // Key the phased VCFs by cohort and chromosome
    ch_phased_vcfs_keyed = ch_phased_vcf_for_imputation
         .map { meta, vcf, index -> ["${meta.chr}_${meta.cohort}", meta, vcf, index] }
         .view { "Phased VCF Keyed: Key ${it[0]}, ID ${it[1].id}" }

    // Enhance the phased_vcf channel to make it more accessible for joining
    ch_phased_vcf_by_cohort = ch_phased_vcf_for_imputation
        .map { meta, vcf, index -> [meta.cohort, meta, vcf, index] }
        .view { "Phased VCF By Cohort: Cohort ${it[0]}, Chr ${it[1].chr}" }

    ch_final_imputed_vcf = Channel.empty()

    // --- Imputation Step (Conditional) ---
    if (params.imputer == 'minimac4') {
        // 1. Create Minimac4 reference panels (.m3vcf)
        CONVERT_REFP(ch_phased_vcf_for_imputation
            .map { meta, vcf, index ->
                if (!meta || !vcf) {
                    error "Missing meta or VCF file for reference panel conversion: meta=${meta}, vcf=${vcf}"
                }
                return [meta, vcf, index]
            }
        )
        
        // Store reference panels in a channel format we can join with
        ch_ref_panels = CONVERT_REFP.out.reference
            .map { meta, msav -> 
                def key = "${meta.chr}_${meta.cohort}"
                log.info "Created reference panel: key=${key}, meta=${meta.id}"
                return [key, meta, msav]
            }
            .view { key, meta, panel -> "Minimac4 Ref Panel: Key ${key}, Meta ${meta.id}, File ${panel}" }
        
        // FIXED APPROACH: Create all possible imputation jobs explicitly
        // Create a list of all VCFs and a list of all reference panels, then create all valid combinations
        // Step 1: Collect all phased VCFs into a list with their metadata and chromosome info
        ch_all_vcfs = ch_phased_vcf_for_imputation
            .map { meta, vcf, index -> 
                [meta.cohort, meta.chr, meta, vcf, index]
            }
            .view { src_cohort, chr, meta, vcf, idx -> 
                "Phased VCF for jobs: Cohort=${src_cohort}, Chr=${chr}" 
            }
            .toList()
        
        // Step 2: Collect all reference panels into a list
        ch_all_ref_panels = ch_ref_panels
            .map { key, meta, panel -> 
                def parts = key.split('_')
                def chr = parts[0]
                def target_cohort = parts[1]
                [target_cohort, chr, meta, panel]
            }
            .view { target_cohort, chr, meta, panel -> 
                "Ref panel for jobs: Target=${target_cohort}, Chr=${chr}" 
            }
            .toList()
        
        // Step 3: Create all valid imputation jobs by cartesian product and filtering
        ch_impute_input_minimac4 = ch_all_vcfs
            .combine(ch_all_ref_panels)
            .map { vcf_list, panel_list ->
                def jobs = []
                
                // For each phased VCF
                vcf_list.each { src_cohort, src_chr, src_meta, src_vcf, src_idx ->
                    // For each reference panel
                    panel_list.each { target_cohort, target_chr, target_meta, target_panel ->
                        // Only create jobs where:
                        // 1. Source and target cohorts are different
                        // 2. Chromosomes match
                        if (src_cohort != target_cohort && src_chr == target_chr) {
                            // Create metadata for the imputation job
                            def meta_impute = src_meta + [target_cohort: target_cohort]
                            
                            log.info "Creating imputation job: ${meta_impute.id} (target: ${target_cohort}), Chr=${src_chr}"
                            
                            // Add job to the list
                            jobs << [meta_impute, src_vcf, src_idx, target_panel]
                        }
                    }
                }
                
                return jobs
            }
            .flatten().buffer(size: 4) // Convert back to a channel of individual jobs
            .view { meta, vcf, idx, panel -> 
                "Imputation job ready: ${meta.id} -> ${meta.target_cohort}" 
            }

        // Run Minimac4 Imputation with the matched inputs
        IMPUTE_MINIMAC4(ch_impute_input_minimac4)
        ch_final_imputed_vcf = IMPUTE_MINIMAC4.out.imputed_vcf

    } else if (params.imputer == 'beagle') {
        // Prepare input for IMPUTE_BEAGLE
        ch_impute_input_beagle = ch_imputation_pairs
            .map { source, target -> ["${source}", source, target] } // Key by source
            .join(ch_phased_vcf_for_imputation.map { m,v,i -> [m.cohort, m, v, i] }, by: 0) // Join source with its phased VCF
            .map { _source_key, _source, target, meta_s, vcf_s, index_s ->
                 def meta_impute = meta_s + [target_cohort: target] // Add target info
                 def ref_key = "${meta_s.chr}_${target}"          // Key to join with target's phased VCF
                 return [ ref_key, meta_impute, vcf_s, index_s ]
            }
            .join(ch_phased_vcfs_keyed, by: 0) // Join with target's *phased VCF* using key [chr_target]
            .map { _ref_key, meta_i, vcf_s, index_s, meta_r, vcf_r, _index_r ->
                 if (meta_i.chr != meta_r.chr) { error "Chromosome mismatch during Beagle input join!" }
                 return [ meta_i, vcf_s, index_s, vcf_r ]
            }

        // Run Beagle Imputation
        IMPUTE_BEAGLE(ch_impute_input_beagle)
        ch_final_imputed_vcf = IMPUTE_BEAGLE.out.imputed_vcf

    } else {
        error "Invalid imputer specified: ${params.imputer}. Choose 'minimac4' or 'beagle'."
    }

    emit:
    imputed_vcf = ch_final_imputed_vcf // Output [meta, vcf, index]
}
