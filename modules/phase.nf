// Process 3: PHASE
// Phase haplotypes using Eagle
process PHASE {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.cohort}/phased", mode: 'copy'
    
    input:
    tuple val(meta), path(bcf), path(bcf_index)

    output:
    tuple val(meta), path("${meta.id}_phased.vcf.gz"), path("${meta.id}_phased.vcf.gz.tbi"), emit: phased_vcf

    script:
    def prefix = "${meta.id}_phased"
    """
    echo "Phasing ${bcf} with Eagle"
    
    # Run Eagle for phasing
    ${params.eagle} \
      --vcf ${bcf} \
      --Kpbwt=10000 \
      --pbwtIters=${params.rounds} \
      --geneticMapFile=${params.gmap} \
      --vcfOutFormat=z \
      --outPrefix ${prefix} \
      --numThreads ${task.cpus}
    
    # Index the phased output
    ${params.tabix} -p vcf ${prefix}.vcf.gz
    
    echo "Phasing completed for ${meta.id}"
    """
} 