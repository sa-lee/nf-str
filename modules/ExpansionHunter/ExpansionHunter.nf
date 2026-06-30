include { call } from './call.nf'
include { kmer_filter } from './kmer_filter.nf'
include { sort_bamlet } from './sort_bamlet.nf'

workflow run_expansion_hunter {
    take:
        sam_bam_ch
    main:
        eh5_results = sam_bam_ch |
            call |
            sort_bamlet 
            //| kmer_filter 
    emit:
        eh5_results  
}