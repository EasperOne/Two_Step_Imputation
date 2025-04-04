// Process 6: PREPARE_R_INPUT
// Extract info files from imputed VCFs for R analysis
process PREPARE_R_INPUT {
    tag "${meta.id} (target: ${meta.target_cohort})"
    publishDir "${params.outdir}/${meta.cohort}/r_input/${meta.target_cohort}", mode: 'copy'
    
    conda "conda-forge::perl=5.26.2"

    input:
    tuple val(meta), path(imputed_vcf)

    output:
    tuple val(meta), path("${meta.id}_imputed_${meta.target_cohort}.all.info.txt"), emit: info_files

    script:
    def output_prefix = "${meta.id}_imputed_${meta.target_cohort}"
    """
    echo "Extracting info from ${imputed_vcf}"
    
    # Run the Perl script to extract info from imputed VCF
    perl ${params.extract_info_pl} ./ ${output_prefix}
    
    echo "Info extraction completed for ${meta.id}"
    """
} 