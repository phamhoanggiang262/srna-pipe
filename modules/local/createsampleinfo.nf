
process CREATESAMPLEINFO {
    tag '$bam'
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

input:
	path samplesheet
	
	output:
	path 'sampleInfo.csv'			,emit: sampleInfo
	
	script:
	"""
	createSampleInfo.R --samplesheet $samplesheet --path $baseDir/results/samtools/sort
	"""	
}

