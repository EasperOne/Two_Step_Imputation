#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Include modules
include { PREPARE_R_INPUT } from '../modules/prepare_r_input'
include { RUN_R_SELECT } from '../modules/run_r_select'
include { FILTER_VCF } from '../modules/filter_vcf'
include { MERGE_VCF } from '../modules/merge_vcf'

workflow POST_IMPUTATION {
    take:
    ch_imputed_vcf // Channel of [meta, vcf] from imputation

    main:
    // First, ensure each input has a VCF index by creating one if needed
    ch_final_imputed_vcf = ch_imputed_vcf
        
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
    
    // Create empty placeholder file if it doesn't exist
    file("${workflow.projectDir}/EMPTY_FILE").text = ""

    // --- Extract keep and filterout lists from R output ---
    ch_parsed_manifest = RUN_R_SELECT.out.filter_manifest
        .flatMap { meta_r, manifest_path ->
            def results = []
            manifest_path.readLines().drop(1).each { line ->
                def fields = line.split('\t')
                if (fields.size() >= 4) {
                    def (chr, type, ref_cohort, file_path) = fields
                    def key = "${chr}_${ref_cohort ?: 'NONE'}"
                    results << [key: key, chr: chr, type: type, ref_cohort: ref_cohort, file: file(file_path)]
                }
            }
            return results
        }
        .view { "Manifest Entry: $it" }

    // Group the manifest entries by key and build filter map
    ch_filter_map = ch_parsed_manifest
        .groupTuple(by: 'key')
        .map { entry ->
            def key = entry.key[0]
            def keep_file = entry.type.contains('keep') ? 
                entry.file[entry.type.indexOf('keep')] : 
                file("${workflow.projectDir}/EMPTY_FILE")
            def filterout_file = entry.type.contains('filterout') ? 
                entry.file[entry.type.indexOf('filterout')] : 
                file("${workflow.projectDir}/EMPTY_FILE")
            [key, filterout_file, keep_file]
        }
        .view { "Filter Files Map: $it" }

    // Prepare imputed VCFs for filtering with proper keys
    ch_imputed_for_filter = ch_imputed_vcf
        .map { meta, vcf, index -> 
            def key = "${meta.chr}_${meta.target_cohort ?: 'NONE'}"
            [key, meta, vcf, index]
        }
        .view { "Imputed for filter: $it" }

    // Join imputed VCFs with filter files
    ch_filter_input = ch_imputed_for_filter
        .join(ch_filter_map, by: 0, remainder: true)
        .map { key, meta, vcf, index, filterout, keep ->
            // Use default empty files if no filter files found
            def filterout_final = filterout ?: file("${workflow.projectDir}/EMPTY_FILE")
            def keep_final = keep ?: file("${workflow.projectDir}/EMPTY_FILE")
            [meta, vcf, index, filterout_final, keep_final]
        }
        .view { "Filter Input: $it" }

    // Run the filter VCF process
    FILTER_VCF(ch_filter_input)
    
    // --- Merge Setup ---
    // Group filtered VCFs by cohort and chromosome for merging
    ch_merge_input = FILTER_VCF.out.filtered_vcf
        .map { meta, vcf, index -> 
            return [ meta.cohort, meta.chr, meta, vcf, index ]
        }
        .groupTuple(by: [0, 1]) // Group by [cohort, chr]
        .map { cohort, chr, metas, vcfs, indices ->
            def meta = [
                id: "${cohort}_chr${chr}_merged",
                cohort: cohort,
                chr: chr
            ]
            return [ meta, vcfs, indices ]
        }
    
    MERGE_VCF(ch_merge_input)

    emit:
    merged_vcf = MERGE_VCF.out.merged_vcf // Final output: [meta, vcf, index]
}
