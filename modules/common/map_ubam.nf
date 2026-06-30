process minimap2_ubam_ont {
    memory { 60.GB * task.attempt }
    cpus 8
    time 96.h
    queue 'long' 
    maxRetries 2
    maxForks 40

    publishDir "output/alignments/", mode: "copy", pattern: "*.{sorted.bam,sorted.bam.bai}"
    
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
    
    tag { sam }
    
    input:
    tuple val(sam), val(type), path(bam)
    
    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam"), path("${sam}_${type}.sorted.bam.bai"), emit:aligned
    tuple val(sam), path(bam), val("SUCCESS"), emit: cleanup_status  // Status for cleanup
    
    script:
    """
    # Stream: BAM -> FASTQ -> minimap2 -> sort -> index
    # No intermediate FASTQ file written to disk
    samtools fastq -T '*' -@ ${task.cpus} ${bam} | \\
        minimap2 -x map-ont -a -t ${task.cpus} --MD ${params.ref_fasta} - | \\
        samtools sort -O bam -@ ${task.cpus} -o ${sam}_${type}.sorted.bam -
    
    # Index the sorted BAM
    samtools index -@ ${task.cpus} ${sam}_${type}.sorted.bam
    """
}

process minimap2_ubam_pacbio {
    memory { 40.GB * task.attempt }
    cpus 8
    maxRetries 2
    maxForks 40

    publishDir "output/alignments/", mode: "copy", pattern: "*.{sorted.bam,sorted.bai}"

    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
    
    tag { sam }
    
    input:
    tuple val(sam), val(type), path(bam)
    
    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam"), path("${sam}_${type}.sorted.bam.bai"), emit: aligned
    tuple val(sam), path(bam), val("SUCCESS"), emit: cleanup_status  // Status for cleanup
    
    script:
    """
    # Stream: BAM -> FASTQ -> minimap2 -> sort -> index
    # No intermediate FASTQ file written to disk
    samtools fastq -T '*' -@ ${task.cpus} ${bam} | \\
        minimap2 -x map-hifi -a -t ${task.cpus} --MD ${params.ref_fasta} - | \\
        samtools sort -O bam -@ ${task.cpus} -o ${sam}_${type}.sorted.bam -
    
    # Index the sorted BAM
    samtools index -@ ${task.cpus} ${sam}_${type}.sorted.bam
    """
}

// Updated Illumina process with increased resources (was already optimized for piping)
process minimap2_ubam_illumina {
    memory { 40.GB * task.attempt }
    cpus 8
    maxRetries 2

    publishDir "output/alignments/", mode: "copy", pattern: "*.{sorted.bam,sorted.bai}"
    
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'

    tag { sam }
    
    input:
    tuple val(sam), val(type), path(bam)
    
    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam"), path("${sam}_${type}.sorted.bam.bai")
    
    script:
    """
    # Stream: BAM -> FASTQ -> minimap2 -> sort -> index
    # No intermediate FASTQ file written to disk
    samtools fastq -T "*" -@ ${task.cpus} ${bam} | \\
        minimap2 -ax sr -t ${task.cpus} ${params.ref_fasta} - | \\
        samtools sort -O bam -@ ${task.cpus} -o ${sam}_${type}.sorted.bam -
    
    # Index the sorted BAM
    samtools index -@ ${task.cpus} ${sam}_${type}.sorted.bam
    """
}