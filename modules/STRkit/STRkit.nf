workflow run_strkit {
    take:
        sam_bam_ch
    main:
        strkit_results = sam_bam_ch |
            strkit
    emit:
        strkit_results  
}

process strkit {
    cpus 1
	memory {'4 GB'}
	time '5 h'
    publishDir "output/strkit/", mode: "copy"

    tag "${sam}_${type}"

    container 'ghcr.io/davidlougheed/strkit:0.24.2'

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}_${type}.vcf")

    script:
    """
    strkit call \
    ${bam} --ref ${params.ref_fasta} --loci ${params.strkit_loci} \
    --vcf ${sam}_${type}.vcf --sample-id ${sam} --seed 123 --no-tsv \
    --hq --realign
    #--processes X --use-hp --incorporate-snvs path/to/dbsnp/00-common_all.vcf.gz 
    """
}

