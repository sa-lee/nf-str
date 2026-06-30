process call {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "output/eh5/", mode: "copy", pattern: "*.{json,vcf}"
    tag { sam }

    container 'quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5'

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}_realigned.bam"), path("${sam}.json"), path("${sam}.vcf")

    script:
    """
        ExpansionHunter \
		--reads ${bam} \
		--reference ${params.ref_fasta} \
		--variant-catalog ${params.eh5_loci} \
		--output-prefix ${sam}
    """
}