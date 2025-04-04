// modules/impute_beagle.nf
// Impute genotypes using Beagle 5.5

process IMPUTE_BEAGLE {
    tag "${meta.id}"
    label 'process_high'
    
    conda "bioconda::beagle=5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/beagle:5.4--hdfd78af_0' :
        'quay.io/biocontainers/beagle:5.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(index)

    output:
    tuple val(meta), path("*_imputed.vcf.gz"), emit: imputed_vcf
    tuple val(meta), path("*_imputed.vcf.gz.tbi"), emit: imputed_index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory.toGiga()
    """
    # Rename the output to avoid collisions
    beagle \\
        gt=${vcf} \\
        out=${prefix}_imputed \\
        nthreads=${task.cpus} \\
        window=${params.window_size} \\
        overlap=${params.overlap_size} \\
        $args

    tabix -p vcf ${prefix}_imputed.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: \$(echo \$(beagle 2>&1 | grep version | sed 's/^.*version //; s/ .*\$//'))
    END_VERSIONS
    """
}