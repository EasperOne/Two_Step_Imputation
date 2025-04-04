#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define default parameter values
params.help = false

// Include modules
include { ALIGNMENT } from './modules/alignment'
include { VCF_TO_BCF } from './modules/vcf_to_bcf'
// Phasing Modules
include { PHASE as PHASE_EAGLE } from './modules/phase' // Rename existing phase module
include { PHASE_SHAPEIT5 } from './modules/phase_shapeit5'
// Conversion Module (if needed after SHAPEIT5)
include { BCF_TO_VCF } from './modules/bcf_to_vcf'
// Reference Panel Module (Only for Minimac4)
include { CONVERT_REFP } from './modules/convert_refp'
// Imputation Modules
include { IMPUTE as IMPUTE_MINIMAC4 } from './modules/impute' // Rename existing impute module
include { IMPUTE_BEAGLE } from './modules/impute_beagle'
// Downstream Modules
include { PREPARE_R_INPUT } from './modules/prepare_r_input'
include { RUN_R_SELECT } from './modules/run_r_select'
include { FILTER_VCF } from './modules/filter_vcf'
include { MERGE_VCF } from './modules/merge_vcf'

// --- Log Pipeline Header & Help ---
def printHeader() {
    log.info """
         Two Step Imputation Pipeline - Nextflow
         ===========================================
         Cohorts CSV   : ${params.cohorts_csv}
         Output Dir    : ${params.outdir}
         Ref FASTA     : ${params.fasta}
         Genetic Map   : ${params.gmap}
         CPUs          : ${params.cpus}
         Minimac Rounds: ${params.rounds}
         ===========================================
         NOTE: This pipeline performs cross-imputation 
         between cohorts. Each cohort is imputed using 
         all other cohorts as reference panels.
         Reference panels are created during the pipeline run.
         ===========================================
         """ .stripIndent()
}

def printHelp() {
    log.info """
    Usage:

    nextflow run <your_repo/main.nf> -profile <standard/test/cluster> [options]

    Required Parameters:
    --cohorts_csv   Path to CSV file describing cohorts (see README).
    --outdir        Directory where results will be saved.
    --fasta         Path to reference genome FASTA file.
    --gmap          Path to genetic map file (for Eagle/SHAPEIT5/Beagle).

    Optional Parameters:
    --phaser        Phasing tool to use ('eagle' or 'shapeit5') (Default: ${params.phaser}).
    --imputer       Imputation tool to use ('minimac4' or 'beagle') (Default: ${params.imputer}).
    --snp_select_r  Path to SNP_Selection.R script (Default: ${params.snp_select_r}).
    --extract_info_pl Path to extract_info_VCF.pl script (Default: ${params.extract_info_pl}).
    --cpus          Base number of CPUs for tasks (Default: ${params.cpus}).
    --rounds        Number of phasing rounds for Minimac4 (Default: ${params.rounds}).
    --chromosomes   Range of chromosomes to process (e.g., \"1..22\") (Default: 1..22).

    Tool Path Parameters (Required if not in PATH/handled by Conda):
    --eagle         Path to Eagle executable.
    --shapeit5      Path to SHAPEIT5 executable.
    --minimac4      Path to Minimac4 executable.
    --beagle        Path to Beagle JAR file.
    --java          Path to Java executable (if not in PATH).
    --bcftools      Path to bcftools executable.
    --tabix         Path to tabix executable.

    Other Options:
    -profile        Configuration profile to use (e.g., standard, test).
    -resume         Resume previous run from cache.
    --help          Show this help message.
    """ .stripIndent()
}

// --- Workflow Definition ---
workflow {

    printHeader()
    
    // Show help message if --help is specified
    if (params.help) {
        printHelp()
        exit 0
    }

    // --- Validate Parameters ---
    // Ensure required files/params exist
    if (!params.fasta || !file(params.fasta).exists()) {
        exit 1, "Reference FASTA file not found: ${params.fasta ?: 'parameter not set'}"
    }
    if (!params.gmap || !file(params.gmap).exists()) {
        exit 1, "Genetic Map file not found: ${params.gmap ?: 'parameter not set'}"
    }
    if (!params.cohorts_csv || !file(params.cohorts_csv).exists()) {
        exit 1, "Cohorts CSV file not found: ${params.cohorts_csv ?: 'parameter not set'}"
    }
    // Optional: Validate other script paths if needed
    // if (!params.snp_select_r || !file(params.snp_select_r).exists()) { exit 1, "..." }
    // if (!params.extract_info_pl || !file(params.extract_info_pl).exists()) { exit 1, "..." }

    // --- Input Channel Creation ---
    // Channel to read Cohorts_Info.csv and find input VCF files per chromosome
    ch_input_vcfs = Channel.fromPath(params.cohorts_csv)
        .splitCsv(header:true, sep:',') // each row is a list of named columns
        .map { row ->
            def cohort = row.Cohort
            def pathway = row.Pathway.trim()
            def prefix = row.Prefix.trim()
            def suffix = row.Suffix.trim()
            // Generate chromosome list (1-22)
            def chromosomes = params.chromosomes ?: (1..22)
            // Create [meta, vcf_file] tuples for each existing VCF
            // collect operator creates a list of tuples for each chromosome in the initiated list
            return chromosomes.collect { chr ->
                def meta = [
                    id: "${cohort}_chr${chr}",
                    cohort: cohort,
                    chr: chr
                ]
                def vcf_path = "${pathway}/${prefix}${chr}${suffix}"
                def vcf_file = file(vcf_path) // Create file object
                if (vcf_file.exists()) {
                    log.info "Found input VCF: ${vcf_file}"
                    return [ meta, vcf_file ]
                } else {
                    return null
                }
            }
        }
        // output nested list of tuples created by the collect operator
        //    [
        //    [ tuple1, tuple2, ..., tuple22 ],  // For cohort1
        //    [ tuple1, tuple2, ..., tuple22 ],  // For cohort2
        //    [ tuple1, tuple2, ..., tuple22 ]   // For cohort3
        //    ]
        .flatMap() // Flatten the list of lists into individual [meta, vcf_path] emissions
        .filter { it != null } // Remove null entries where files didn't exist
        .ifEmpty { exit 1, "No input VCF files found based on ${params.cohorts_csv}. Please check paths and naming convention (e.g., path/prefixCHRsuffix)." }

    // Pipeline steps
    ALIGNMENT(ch_input_vcfs)
    VCF_TO_BCF(ALIGNMENT.out.vcf)

    // --- Phasing Step (Conditional) ---
    if (params.phaser == 'eagle') {
        PHASE_EAGLE(VCF_TO_BCF.out.bcf)
        ch_phased_vcf_for_imputation = PHASE_EAGLE.out.phased_vcf // Eagle outputs VCF.gz
    } else if (params.phaser == 'shapeit5') {
        PHASE_SHAPEIT5(VCF_TO_BCF.out.bcf)
        BCF_TO_VCF(PHASE_SHAPEIT5.out.phased_bcf) // Convert SHAPEIT5 BCF to VCF
        ch_phased_vcf_for_imputation = BCF_TO_VCF.out.converted_vcf
    } else {
        exit 1, "Invalid phaser specified: ${params.phaser}. Choose 'eagle' or 'shapeit5'."
    }
    ch_phased_vcf_for_imputation.view { "Phased VCF ready for Imputation/RefPanel: ${it[0].id}" }

    // --- Imputation Setup (Part 1: Reference Panels) ---
    // Reference Panel generation OR selection depends on imputer choice
    // Create channels for imputation pairs and keyed phased VCFs (needed by both imputers)
    ch_cohort_names = ch_input_vcfs
        .map { meta, _vcf -> meta.cohort }
        .unique()
        .collect()
        .view { "Cohort Names: $it" }

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

    // Key the phased VCFs (now consistently VCF.gz) by cohort and chromosome
    ch_phased_vcfs_keyed = ch_phased_vcf_for_imputation
         .map { meta, vcf, index -> ["${meta.chr}_${meta.cohort}", meta, vcf, index] }
         .view { "Phased VCF Keyed: Key ${it[0]}, ID ${it[1].id}" }

    // Enhance the phased_vcf channel to make it more accessible for joining
    ch_phased_vcf_by_cohort = ch_phased_vcf_for_imputation
        .map { meta, vcf, index -> [meta.cohort, meta, vcf, index] }
        .view { "Phased VCF By Cohort: Cohort ${it[0]}, Chr ${it[1].chr}" }

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
            
        // Modify the approach to avoid using cross operation
        // First, make the imputation pairs have source cohort as key
        ch_impute_source_keys = ch_imputation_pairs
            .map { source, target -> [source, [source, target]] }
        
        // Now join with phased VCFs (which are already keyed by cohort)
        ch_impute_source_joined = ch_impute_source_keys
            .join(ch_phased_vcf_by_cohort)
            .map { source, pair_info, meta_s, vcf_s, idx_s ->
                def target = pair_info[1]
                def meta_impute = meta_s + [target_cohort: target]
                def ref_key = "${meta_s.chr}_${target}" // Key to find reference panel
                
                log.info "Created imputation job: Source=${source}, Target=${target}, Chr=${meta_s.chr}, RefKey=${ref_key}"
                
                // Return [ref_key, meta with target info, vcf, index]
                return [ref_key, meta_impute, vcf_s, idx_s]
            }
            
        // Simplify the channel approach - remove the tap and count operations
        ch_impute_input_minimac4 = ch_impute_source_joined
            .join(ch_ref_panels, by: 0) // join on ref_key
            .map { ref_key, meta_impute, vcf_s, idx_s, meta_r, ref_panel ->
                log.info "Creating imputation job: ${meta_impute.id} with ref panel from ${meta_r.id}"
                return [meta_impute, vcf_s, idx_s, ref_panel]
            }

        // 3. Run Minimac4 Imputation with the matched inputs
        IMPUTE_MINIMAC4(ch_impute_input_minimac4)
        ch_final_imputed_vcf = IMPUTE_MINIMAC4.out.imputed_vcf

    } else if (params.imputer == 'beagle') {
        // 1. Prepare input for IMPUTE_BEAGLE
        //    Join [Source, Target] pairs with phased VCF of Source
        //    Then join with *phased VCF* of Target (used as ref panel)
        ch_impute_input_beagle = ch_imputation_pairs
            .map { source, target -> ["${source}", source, target] } // Key by source
            .join( ch_phased_vcf_for_imputation.map { m,v,i -> [m.cohort, m, v, i] }, by: 0 ) // Join source with its phased VCF
            .map { _source_key, _source, target, meta_s, vcf_s, index_s ->
                 def meta_impute = meta_s + [target_cohort: target] // Add target info
                 def ref_key = "${meta_s.chr}_${target}"          // Key to join with target's phased VCF
                 return [ ref_key, meta_impute, vcf_s, index_s ]
            }
            .join(ch_phased_vcfs_keyed, by: 0) // Join with target's *phased VCF* using key [chr_target]
            .map { _ref_key, meta_i, vcf_s, index_s, meta_r, vcf_r, _index_r ->
                 // Final input: [meta_for_imputation, source_vcf, source_index, ref_panel_vcf]
                 // Ensure meta_r matches expectations if needed (same chr etc)
                 if (meta_i.chr != meta_r.chr) { error "Chromosome mismatch during Beagle input join!" }
                 return [ meta_i, vcf_s, index_s, vcf_r ]
            }

        // 2. Run Beagle Imputation
        IMPUTE_BEAGLE(ch_impute_input_beagle)
        ch_final_imputed_vcf = IMPUTE_BEAGLE.out.imputed_vcf

    } else {
        exit 1, "Invalid imputer specified: ${params.imputer}. Choose 'minimac4' or 'beagle'."
    }

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

    // Prepare imputed VCFs for joining (Input is ch_final_imputed_vcf)
    ch_imputed_vcfs_keyed = ch_final_imputed_vcf
        .filter { meta, _vcf, _index -> meta?.chr != null && meta?.target_cohort != null }
        .map { meta, vcf, _index -> // Handle VCF + Index tuple
            def key_keep = "${meta.chr}_${meta.target_cohort}"
            def key_filterout = meta.chr
            return [ key_keep, key_filterout, meta, vcf ]
        }
        .view { tuple -> 
            // Unpack the tuple explicitly
            def (kk, kf, m, v) = tuple
            return "Imputed VCF Keyed: KeepKey ${kk}, FilterKey ${kf}, ID ${m.id}"
        }

    // Join imputed VCFs with keep lists and then with filterout lists
    ch_filter_input = ch_imputed_vcfs_keyed
        .join(ch_keep_lists, by: 0) // Join by key_keep
        .view { tuple ->
            // Unpack the tuple explicitly
            def (kk, kf, m, v, kl) = tuple
            return "Joined Keep: KeepKey ${kk}, FilterKey ${kf}, ID ${m.id}, KeepFile ${kl}"
        }
        .map { _key_keep, key_filterout, meta, vcf, keep_list ->
            // Re-key for joining with filterout list
            return [ key_filterout, meta, vcf, keep_list ]
        }
        .join(ch_filterout_lists, by: 0) // Join by key_filterout
        .view { kf, m, v, kl, fl -> "Joined Filterout: FilterKey ${kf}, ID ${m.id}, KeepFile ${kl}, FilterFile ${fl}" }
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

    workflow.onComplete {
        log.info ( "Pipeline Complete" )
    }

    workflow.onError {
        log.error "Pipeline Failed!"
        if (workflow.errorMessage) {
            log.error "Error message: ${workflow.errorMessage}"
        }
        if (workflow.exitStatus) {
            log.error "Exit status: ${workflow.exitStatus}"
        }
    }
}