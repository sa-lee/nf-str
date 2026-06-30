#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.id = 'v3'
params.ref_fasta = '/stornext/Bioinf/data/lab_bahlo/ref_db/human/hg38/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa'
params.illumina_ref_fasta = '/stornext/Bioinf/data/lab_bahlo/ref_db/human/hg38/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa'

// Manifest for downloading and aligning data (if not skipping download/align steps)
params.download_manifest = '/vast/scratch/users/reid.j/nf-str-run/v3/long_submission.tsv'
params.cleanup_intermediates = false

// Alternate manifest for already downloaded/aligned data (if skipping download/align steps)
params.aligned_manifest = '/vast/projects/reidj-project/nf-str/test_manifest.tsv'
params.skip_download_align = true

// Repeat catalogues
params.atarva_loci = "${projectDir}/repeat_catalogs/STRchive-disease-loci-hg38_atarva.bed.gz"
params.eh5_loci = "${projectDir}/repeat_catalogs/STRchive-disease-loci-hg38_expansionhunter.json"
params.longtr_loci = "${projectDir}/repeat_catalogs/STRchive-disease-loci-hg38_longtr.bed"
params.straglr_loci = "${projectDir}/repeat_catalogs/STRchive-disease-loci-hg38_straglr.bed"
params.strkit_loci = "${projectDir}/repeat_catalogs/STRchive-disease-loci-hg38_strkit.bed"
params.trgt_loci = "${projectDir}/repeat_catalogs/STRchive-disease-loci-hg38_trgt.bed"

include { path; read_tsv; date_ymd } from './modules/functions'
include { index_bam } from './modules/common/sort_bam.nf'
include { download_s3_files } from './modules/common/download_s3_files.nf'
include { minimap2_ubam_illumina; minimap2_ubam_ont; minimap2_ubam_pacbio } from './modules/common/map_ubam.nf'

//ExpansionHunter
include { run_expansion_hunter } from './modules/ExpansionHunter/ExpansionHunter.nf'

//ExpansionHunterDeNovo
include { run_expansion_hunter_denovo } from './modules/ExpansionHunterDenovo/ExpansionHunterDenovo.nf'

//Scatter
include { run_scattr } from './modules/ScatTR/ScatTR.nf'

//Straglr
include { run_straglr as run_straglr } from './modules/Straglr/Straglr.nf'

//LongTR
include { run_longtr } from './modules/LongTR/LongTR.nf'

//STRkit
include { run_strkit } from './modules/STRkit/STRkit.nf'

// atarva
include { run_atarva as run_atarva } from './modules/atarva/atarva.nf'

// TRGT
include { run_trgt } from './modules/TRGT/TRGT.nf'

workflow {
    if (params.skip_download_align && params.aligned_manifest) {
        // Use alternate manifest for already downloaded/aligned data
        aligned_manifest = read_tsv(path(params.aligned_manifest), ['sample', 'type', 'bam_path'])
        
        all_samples_with_index = Channel
            .from(aligned_manifest)
            .map { record -> 
                def sample = record.sample
                def type = record.type
                def bam_file = file(record.bam_path)

                // Auto-detect index file in same directory
                def index_file
                if (bam_file.toString().endsWith('.cram')) {
                    index_file = file(bam_file.toString() + '.crai')
                } else {
                    index_file = file(bam_file.toString() + '.bai')
                }
                
                tuple(sample, type, bam_file, index_file)
            }          
    } else {
        // Original workflow logic (download and align)
        manifest = read_tsv(path(params.manifest), ['sample', 'type', 'url', 'align'])

        s3_input_ch = Channel
            .from(manifest)
            .map { record -> 
                def sample = record.sample
                def type = record.type
                def bam_url = record.url 
                def align = record.align?.toLowerCase()?.trim() == 'yes'
                tuple(sample, type, bam_url, align)
            }

        downloaded = download_s3_files(s3_input_ch)
        
        alignment_check = downloaded
            .branch { sample, type, bam_file, needs_align ->
                to_align: needs_align
                    return tuple(sample, type, bam_file)
                already_aligned: !needs_align
                    return tuple(sample, type, bam_file)
            }
        
        unaligned_by_type = alignment_check.to_align
            .branch { sample, type, bam_file ->
                illumina: type == 'illumina'
                    return tuple(sample, type, bam_file)
                ont: type == 'ont'
                    return tuple(sample, type, bam_file)
                pacbio: type == 'pacbio'
                    return tuple(sample, type, bam_file)
            }
      
        illumina_aligned = minimap2_ubam_illumina(unaligned_by_type.illumina)
        ont_aligned = minimap2_ubam_ont(unaligned_by_type.ont)
        pacbio_aligned = minimap2_ubam_pacbio(unaligned_by_type.pacbio)

        ont_aligned.cleanup_status
            .filter { sample_id, bam, status -> status == "SUCCESS" }
            | conditional_cleanup_ont

        pacbio_aligned.cleanup_status
            .filter { sample_id, bam, status -> status == "SUCCESS" }
            | conditional_cleanup_pacbio  
        
        all_aligned = illumina_aligned
            .mix(ont_aligned.aligned, pacbio_aligned.aligned)
        
        already_aligned_samples = alignment_check.already_aligned
        
        already_aligned_samples
            .map { sample, type, bam_file ->
                def index_file
                def index_exists = false
                
                if (bam_file.toString().endsWith('.cram')) {
                    index_file = file(bam_file.toString() + '.crai')
                } else {
                    index_file = file(bam_file.toString() + '.bai')
                }
                
                if (index_file.exists()) {
                    index_exists = true
                }
                
                tuple(sample, type, bam_file, index_file, index_exists)
            }
            .branch { sample, type, bam, index, exists ->
                has_index: 
                    exists == true
                    return tuple(sample, type, bam, index)
                needs_index: 
                    exists == false
                    return tuple(sample, type, bam)
            }
            .set { indexed_check }
        
        generated_indexes = index_bam(indexed_check.needs_index)
        
        all_samples_with_index = all_aligned
            .mix(indexed_check.has_index, generated_indexes)
    }
    // Common downstream processing (same for both entry points)
    all_samples_with_index
        .branch { sample, type, bam, index ->
            illumina: 
                type == 'illumina'
                return tuple(sample, type, bam, index)
            ont: 
                type == 'ont'
                return tuple(sample, type, bam, index)
            pacbio: 
                type == 'pacbio'
                return tuple(sample, type, bam, index)
        }
        .set { samples }
    
    illumina_results = run_illumina(samples.illumina)
    ont_results = run_ont(samples.ont)
    pacbio_results = run_pacbio(samples.pacbio)
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

workflow run_illumina {
    take:
        sample_ch
    main:
        eh_results = sample_ch  | run_expansion_hunter
        //ehdn_results = sample_ch | run_expansion_hunter_denovo
        //scattr_results = sample_ch | run_scattr
            
    emit:
        eh_results
        //ehdn_results
        //scattr_results
}

workflow run_ont {
    take:
        sample_ch
    main:
        atarva_results = sample_ch | run_atarva
        straglr_results = sample_ch | run_straglr
        longtr_results = sample_ch | run_longtr
        strkit_results = sample_ch | run_strkit
        
    emit:
        straglr_results
        longtr_results
        atarva_results
        strkit_results
}

workflow run_pacbio {
    take:
        sample_ch
    main:
        atarva_results = sample_ch | run_atarva
        trgt_results = sample_ch | run_trgt
        straglr_results = sample_ch | run_straglr
        longtr_results = sample_ch | run_longtr
        strkit_results = sample_ch | run_strkit
        
    emit:
        straglr_results
        longtr_results
        atarva_results
        strkit_results
}



process conditional_cleanup_ont {
    when params.cleanup_intermediates == true
    
    input:
    tuple val(sam), path(bam), val(status)
    
    script:
    """
    echo "Symlink: ${bam}"
    echo "Realpath: \$(realpath ${bam})"
    #rm -f \$(realpath ${bam})
    truncate -s 0 \$(realpath ${bam})
    """
}


process conditional_cleanup_pacbio {
    when params.cleanup_intermediates == true
    
    input:
    tuple val(sam), path(bam), val(status)
    
    script:
    """
    echo "Symlink: ${bam}"
    echo "Realpath: \$(realpath ${bam})"
    truncate -s 0 \$(realpath ${bam})
    """
}
