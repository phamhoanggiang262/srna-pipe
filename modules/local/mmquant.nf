process MMQUANT {
    tag "$meta.id"
    label 'process_long'
    publishDir  "$baseDir/results/mmquant", mode: 'copy'
    
    conda "bioconda::mmquant=1.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmquant:1.0.5--hd03093a_0':
        'quay.io/biocontainers/mmquant:1.0.5--hd03093a_0' }"
       
    input:
    tuple val(meta), path (bam), path (annotation)
	
    output:
    tuple val(meta), path ('count_data.tsv')    				, emit: count_matrix
    tuple val(meta), path ('count_report.tsv')					, emit: count_report
    path "versions.yml"                         				, emit: versions
	
    when:
    task.ext.when == null || task.ext.when
	
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    mmquant \\
      -a $annotation \\
      -r $bam \\
      -o count_data.tsv \\
      -O count_report.tsv \\
       $args
	
	
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmquant: \$(mmquant -v 2>&1 | sed -e "s/mmquant version //g")
    END_VERSIONS
    """
}




