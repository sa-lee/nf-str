workflow run_trgt {
    take:
        sam_bam_ch
    main:
        trgt_results = sam_bam_ch |
            trgt
    emit:
        trgt_results  
}

process trgt {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "output/trgt/", mode: "copy", saveAs: { filename ->  filename.replaceAll("${sam}\\.", "${sam}_${type}.")}


    tag "${sam}_${type}"

    container 'quay.io/biocontainers/trgt:4.1.0--h9ee0642_0'

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}.vcf.gz"), path("${sam}.spanning.bam")

    script:
    """
    trgt genotype --genome ${params.ref_fasta} \
                  --repeats ${params.trgt_loci} \
                  --reads ${bam} \
                  --output-prefix ${sam}
    """
}

'''
./trgt plot --genome example/reference.fasta \
       --repeats example/repeat.bed \
       --vcf sample.sorted.vcf.gz \
       --spanning-reads sample.spanning.sorted.bam \
       --repeat-id TR1 \
       --image TR1.svg

trgt plot -g /stornext/Bioinf/data/lab_bahlo/ref_db/human/hg38/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa -b /vast/projects/reidj-project/nf-str/catalogues/STRchive-disease-loci.hg38.TRGT.bed -i SCA27B_FGF14  --vcf /vast/projects/bahlo_longstr/nf-str-output/trgt/NA21110_pacbio.vcf.gz -r /vast/scratch/users/reid.j/trgt/NA21110_pacbio.spanning.sorted.bam
'''