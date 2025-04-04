#!/usr/bin/env nextflow

// Include modules
include { PREPARE_R_INPUT } from '../modules/prepare_r_input'
include { RUN_R_SELECT } from '../modules/run_r_select'
include { FILTER_VCF } from '../modules/filter_vcf'
include { MERGE_VCF } from '../modules/merge_vcf'

workflow POST_IMPUTATION {
    take:
    ch_final_imputed_vcf // Channel of [meta, vcf, index] from imputation

    main:
    // --- Downstream Analysis --- (Starts from ch_final_imputed_vcf)
    PREPARE_R_INPUT(ch_final_imputed_vcf)
    ch_r_input = PREPARE_R_INPUT.out.info_files
        .map { meta, info_file -> [ meta.cohort, meta, info_file ] }
        .groupTuple(by: 0)
        .map { cohort, _meta_list, info_files_list ->
            def meta_r = [ id: "${cohort}_R_analysis", cohort: cohort ]
            return [ meta_r, info_files_list ]
        }
        .view { m, files -> "R Input Ready: ${m.id} with ${files.size()} info files" }

    RUN_R_SELECT(ch_r_input)

    // --- Filter Setup ---
    // Process the manifest file from RUN_R_SELECT
    ch_parsed_manifest = RUN_R_SELECT.out.filter_manifest
        .flatMap { meta_r, manifest_path ->
            manifest_path.readLines().drop(1) // Read lines, skip header
                .collect { line ->
                    def (chr, type, ref_cohort, file_path) = line.split('\t')
                    // Create keys for joining
                    def key_keep = "${chr}_${ref_cohort}" // e.g., "10_Cohort_test_2"
                    def key_filterout = chr              // e.g., "10"
                    // Return structured data including original meta from R run
                    return [ meta_r: meta_r, chr: chr, type: type, ref_cohort: ref_cohort, file: file(file_path), key_keep: key_keep, key_filterout: key_filterout ]
                }
        }
        .view { "Parsed Manifest: Chr ${it.chr}, Type ${it.type}, Ref ${it.ref_cohort ?: '-'}, Path ${it.file}" }

    // Separate keep and filterout lists
    ch_keep_lists = ch_parsed_manifest
        .filter { it.type == 'keep' }
        .map { [ it.key_keep, it.file ] } // Key: "${chr}_${ref_cohort}"
        .unique { it[0] } // Ensure unique key per keep list
        .view { "Keep List Channel: Key ${it[0]}, File ${it[1]}" }

    ch_filterout_lists = ch_parsed_manifest
        .filter { it.type == 'filterout' }
        .map { [ it.key_filterout, it.file ] } // Key: "${chr}"
        .unique { it[0] } // Ensure unique key per filterout list (one per chromosome)
        .view { "Filterout List Channel: Key ${it[0]}, File ${it[1]}" }

    // Prepare imputed VCFs for joining
    ch_imputed_vcfs_keyed = ch_final_imputed_vcf
        .filter { meta, _vcf, _index -> meta?.chr != null && meta?.target_cohort != null }
        .map { meta, vcf, _index -> // Handle VCF + Index tuple
            def key_keep = "${meta.chr}_${meta.target_cohort}"
            def key_filterout = meta.chr
            return [ key_keep, key_filterout, meta, vcf ]
        }
        .view { tuple -> 
            // Unpack the tuple explicitly
            def (kk, kf, m, _v) = tuple
            return "Imputed VCF Keyed: KeepKey ${kk}, FilterKey ${kf}, ID ${m.id}"
        }

    // Join imputed VCFs with keep lists and then with filterout lists
    ch_filter_input = ch_imputed_vcfs_keyed
        .join(ch_keep_lists, by: 0) // Join by key_keep
        .view { tuple ->
            // Unpack the tuple explicitly
            def (kk, kf, m, _v, kl) = tuple
            return "Joined Keep: KeepKey ${kk}, FilterKey ${kf}, ID ${m.id}, KeepFile ${kl}"
        }
        .map { _key_keep, key_filterout, meta, vcf, keep_list ->
            // Re-key for joining with filterout list
            return [ key_filterout, meta, vcf, keep_list ]
        }
        .join(ch_filterout_lists, by: 0) // Join by key_filterout
        .view { kf, m, _v, kl, fl -> "Joined Filterout: FilterKey ${kf}, ID ${m.id}, KeepFile ${kl}, FilterFile ${fl}" }
        .map { _key_filterout, meta, vcf, keep_list, filterout_list ->
            // Final mapping to FILTER_VCF input format
            return [ meta, vcf, filterout_list, keep_list ]
        }

    FILTER_VCF(ch_filter_input)
    
    // --- Merge Setup ---
    // Group filtered VCFs by cohort and chromosome for merging
    ch_merge_input = FILTER_VCF.out.filtered_vcf
        .filter { meta, _vcf, _index -> meta != null && meta.cohort != null && meta.chr != null }
        .map { meta, vcf, index -> 
            return [ meta.cohort, meta.chr, meta, vcf, index ]
        }
        .groupTuple(by: [0, 1]) // Group by [cohort, chr]
        .map { _cohort, _chr, _metas, vcfs, indices ->
            def meta = [
                id: "${_cohort}_chr${_chr}_merged",
                cohort: _cohort,
                chr: _chr
            ]
            return [ meta, vcfs, indices ]
        }
    
    MERGE_VCF(ch_merge_input)

    emit:
    merged_vcf = MERGE_VCF.out.merged_vcf // Final output: [meta, vcf, index]
}
