// Process 5: IMPUTE
// Impute genotypes using minimac4 with reference panels created from other cohorts
process IMPUTE {
    tag "${meta.id} (target: ${meta.target_cohort})"
    publishDir "${params.outdir}/${meta.cohort}/imputed/${meta.target_cohort}", mode: 'copy'

    input:
    tuple val(meta), path(phased_vcf), path(phased_vcf_index), path(ref_panel)

    output:
    tuple val(meta), path("${meta.id}_imputed_${meta.target_cohort}.vcf.gz"), emit: imputed_vcf

    script:
    // Define output prefix here to use it in the output block
    output_prefix = "${meta.id}_imputed_${meta.target_cohort}"
    """
    echo "Starting imputation for ${meta.id} (Target: ${meta.target_cohort}) using ${task.cpus} CPUs"
    echo "Reference panel file: ${ref_panel}"
    echo "Input VCF: ${phased_vcf}"

    echo "Minimac4 Version:"
    ${params.minimac4} --version
    echo "Tabix Version:"
    ${params.tabix} --version || true

    # Run minimac4 imputation
    ${params.minimac4} \
      ${ref_panel} \
      ${phased_vcf} \
      --output ${output_prefix}.vcf.gz \
      -f GT \
      --all-typed-sites \
      --threads ${task.cpus}
    
    echo "Imputation completed for ${meta.id} with reference panel from ${meta.target_cohort}"
    """
} 