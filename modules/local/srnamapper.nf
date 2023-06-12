process SRNAMAPPER {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::srnamapper=1.0.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/srnamapper:1.0.8--h7132678_1':
        'biocontainers/srnamapper:1.0.8--h7132678_1' }"

        
    input:
    tuple val(meta) , path(reads), path(index)

    output:
    tuple val(meta), path("*.sam")		, emit: sam
    path "versions.yml"           		, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    gunzip -c $reads > ${prefix}.fastq
    
    srnaMapper \\
    	-r ${prefix}.fastq \\
    	-g \$INDEX \\
    	-o ${prefix}.sam \\
    	$args 
    	 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        srnamapper: \$(srnaMapper -v | sed -e "s/srnaMapper v//g") 
    END_VERSIONS
    """
}
