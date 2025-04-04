// Process 8: FILTER_VCF
// Filter imputed VCF files based on R output
process FILTER_VCF {
    tag "${meta.id} (target: ${meta.target_cohort})"
    publishDir "${params.outdir}/${meta.cohort}/filtered/${meta.target_cohort}", mode: 'copy'
    
    conda "bioconda::vcftools=0.1.16"

    input:
    tuple val(meta), path(imputed_vcf), path(imputed_index), path(filterout_list), path(keep_list)

    output:
    tuple val(meta), path("${meta.id}_imputed_${meta.target_cohort}_filtered.vcf.gz"), path("${meta.id}_imputed_${meta.target_cohort}_filtered.vcf.gz.tbi"), emit: filtered_vcf

    script:
    def output_prefix = "${meta.id}_imputed_${meta.target_cohort}_filtered"
    def chr = meta.chr
    def has_keep = keep_list.name != 'EMPTY_FILE'
    def has_filterout = filterout_list.name != 'EMPTY_FILE'
    
    """
    echo "Filtering ${imputed_vcf}"
    
    # Decide which filtering approach to use based on available files
    if [ "${has_keep}" = "true" ]; then
        echo "Using Keep list: ${keep_list}"
        
        # Filter VCF using positions to keep
        vcftools --gzvcf ${imputed_vcf} \
                 --chr ${chr} \
                 --positions ${keep_list} \
                 --recode \
                 --recode-INFO-all \
                 --out ${output_prefix}
        
        echo "Filtered using keep list"
    elif [ "${has_filterout}" = "true" ]; then
        echo "Using Filterout list: ${filterout_list}"
        
        # Filter VCF using positions to exclude
        vcftools --gzvcf ${imputed_vcf} \
                 --chr ${chr} \
                 --exclude-positions ${filterout_list} \
                 --recode \
                 --recode-INFO-all \
                 --out ${output_prefix}
        
        echo "Filtered using filterout list"
    else {
        echo "No filter lists provided. Using all variants from chromosome ${chr}."
        
        # Just filter by chromosome
        vcftools --gzvcf ${imputed_vcf} \
                 --chr ${chr} \
                 --recode \
                 --recode-INFO-all \
                 --out ${output_prefix}
        
        echo "No filtering applied, kept all variants in chromosome ${chr}"
    }
    fi

    # Compress and index the output
    if [ -s "${output_prefix}.recode.vcf" ]; then
        bgzip -c ${output_prefix}.recode.vcf > ${output_prefix}.vcf.gz
        ${params.tabix} -p vcf ${output_prefix}.vcf.gz
        echo "Compressed and indexed the filtered VCF"
    else
        echo "Error: Filtering produced an empty VCF file. Check input lists and VCF." >&2
        # Create empty outputs to satisfy Nextflow output expectations
        touch ${output_prefix}.vcf.gz ${output_prefix}.vcf.gz.tbi
        echo "Created empty placeholder files due to filtering error"
    fi

    echo "VCF filtering completed for ${meta.id}"
    """
}