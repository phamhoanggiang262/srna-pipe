process MMANNOT {
    tag "$meta.id"
    label 'process_medium'
    publishDir  "$baseDir/results/mmannot", mode: 'copy'

    conda "bioconda::mmannot=1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmannot:1.1--hd03093a_0 ':
        'quay.io/biocontainers/mmannot:1.1--hd03093a_0 ' }"

    input:
    tuple val(meta), path(bam), path (annotation)
    path configfile

    output:
    path 'annotation_file.tsv'				, emit: annotation_report
    path 'statistics.txt'				, emit: out_stats
    path "versions.yml"           			, emit: versions
	
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """ 
    mmannot \\
    	-a $annotation \\
    	-r $bam \\
    	-o annotation_file.tsv \\
    	-c $configfile \\
    	2> statistics.txt
    	
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
   	mmannot : \$(mmannot -v 2>&1 | sed "s/mmannot v//g")
    END_VERSIONS
    
    """
}
