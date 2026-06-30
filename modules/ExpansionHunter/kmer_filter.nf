process kmer_filter {
    
    publishDir "output/eh5/kmer_filter/", mode: "copy"
    tag { sam }


    input:
    tuple val(sam), path(bamlet_srt), path(bamlet_bai), path(vcf)

    output:
    tuple val(sam), path("${sam}_validated.vcf"), path("${sam}_validated.vcf.tbi")
    
    script:
    out_vcf = "${sam}_validated.vcf"
    out_vcf_tbi = "${sam}_validated.vcf.tbi"
    """
    python kmer_filter.py --bam ${bamlet_srt} --vcf ${vcf} --catalog ${params.eh5_loci} --auto --keep_lowdepth \
-o ./${sam}
    tabix -p vcf ${out_vcf}
    """
}