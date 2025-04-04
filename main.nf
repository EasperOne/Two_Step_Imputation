#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define default parameter values
params.help = false

// Include modules (kept for reference, actual usage is in subworkflows)
include { ALIGNMENT } from './modules/alignment'
include { VCF_TO_BCF } from './modules/vcf_to_bcf'
include { PHASE as PHASE_EAGLE } from './modules/phase'
include { PHASE_SHAPEIT5 } from './modules/phase_shapeit5'
include { BCF_TO_VCF } from './modules/bcf_to_vcf'
include { CONVERT_REFP } from './modules/convert_refp'
include { IMPUTE as IMPUTE_MINIMAC4 } from './modules/impute'
include { IMPUTE_BEAGLE } from './modules/impute_beagle'
include { PREPARE_R_INPUT } from './modules/prepare_r_input'
include { RUN_R_SELECT } from './modules/run_r_select'
include { FILTER_VCF } from './modules/filter_vcf'
include { MERGE_VCF } from './modules/merge_vcf'

// Include subworkflows
include { PRE_PHASE } from './subworkflows/pre_phase'
include { PHASE_WORKFLOW } from './subworkflows/phase_workflow'
include { IMPUTATION } from './subworkflows/imputation'
include { POST_IMPUTATION } from './subworkflows/post_imputation'

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

    Other Options:
    -profile        Configuration profile to use (e.g., standard, test).
    -resume         Resume previous run from cache.
    --help          Show this help message.
    """ .stripIndent()
}

// --- Main Workflow Definition ---
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

    // --- Call Subworkflows ---
    
    // 1. Alignment to BCF conversion
    ch_bcf = PRE_PHASE(ch_input_vcfs)
    
    // 2. Phasing (using either eagle or shapeit5)
    ch_phased_vcf = PHASE_WORKFLOW(ch_bcf)
    
    // 3. Imputation
    ch_imputed_vcf = IMPUTATION(ch_phased_vcf, ch_input_vcfs)
    
    // 4. Post-imputation analysis (R processing, filtering, merging)
    POST_IMPUTATION(ch_imputed_vcf)

    // --- Workflow completion handlers ---
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