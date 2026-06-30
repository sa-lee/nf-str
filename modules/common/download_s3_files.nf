process download_s3_files {
    tag "${sample}_${type}"

    publishDir "output/alignments/", mode: 'copy', enabled: "${ !align }"

    module 'awscli'
    container null
    //container 'quay.io/biocontainers/awscli:1.8.3--py35_0'

    cpus 4
    memory '4 GB'
    time '12h'
    maxForks 2  // Limit to N concurrent downloads

    input:
    tuple val(sample), val(type), val(s3_uri), val(align)

    output:
    tuple val(sample), val(type), path("${filename}"), val(align)
    //path("${filename}"), emit: to_delete

    script:
    filename = s3_uri.tokenize('/')[-1]
    """
    # Set AWS CLI config via environment variables (no config file modification)
    export AWS_MAX_CONCURRENT_REQUESTS=50
    export AWS_MAX_BANDWIDTH=500MB/s
    export AWS_MULTIPART_THRESHOLD=64MB
    export AWS_MULTIPART_CHUNKSIZE=16MB
    export AWS_NO_SIGN_REQUEST=YES

    # Download file from S3
    aws s3 cp --no-sign-request ${s3_uri} ${filename}

    # Verify download
    if [ ! -f "${filename}" ]; then
        echo "Error: Failed to download ${filename}" >&2
        exit 1
    fi

    echo "Successfully downloaded ${filename} (\$(du -h ${filename} | cut -f1))"
    """
}
