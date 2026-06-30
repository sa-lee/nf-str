params.caller = 'atarva'
params.atarva_loci = "${projectDir}/catalogues/STRchive-disease-loci.hg38.atarva.sorted.bed.gz"  

workflow run_atarva {
    take:
        sam_bam_ch
    main:
        atarva_results = sam_bam_ch |
            atarva
    emit:
        atarva_results  
}

process atarva {
    cpus 2
    memory {'5 GB'}
    time '5 h'
    publishDir "output/atarva/", mode: "copy", saveAs: { filename ->
        filename.replaceAll("${sam}\\.", "${sam}_${type}.")
    }

    container 'dhaksnamoorthy/atarva:v0.3.1'

    tag "${sam}_${type}"

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}.vcf")

    script:
    def input_format = bam.name.endsWith('.cram') ? 'cram' : 'bam'
    """
    mkdir -p ./cache
    export REF_PATH=${params.ref_fasta}
    export REF_CACHE='./cache/%2s/%2s/%s'

    atarva --fasta ${params.ref_fasta} \\
           --bam ${bam} \\
           --format ${input_format} \\
           --regions ${params.atarva_loci} \\
           --vcf ${sam}.vcf \\
           --threads ${task.cpus}
    """
}


