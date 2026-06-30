process index_bam {
    cpus 8
    memory '10 GB'
    time '1 hours'
    tag { sam }

    publishDir "output/alignments/", mode: "copy", pattern: "*.{bai,crai}"
    
    input:
    tuple val(sam), val(type), path(bam)

    output:
    tuple val(sam), val(type),  path(bam), path("${bam}*")
    
    script:
    """
    samtools index --threads ${task.cpus} ${bam} 
    """
}

process sort_bam {
    cpus 8
    memory '10 GB'
    time '1 hours'
    tag { sam }
    
    input:
    tuple val(sam), val(type), path(bam)

    output:
    tuple val(sam), val(type), path("*.sorted.bam")
    
    script:
    def bam_sorted = bam.replaceAll('.bam', '.sorted.bam')
    """
    samtools sort --threads ${task.cpus} -o ${bam_sorted} ${bam}
    """
}
