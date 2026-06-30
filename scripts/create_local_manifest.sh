#!/usr/bin/env bash

# Usage: bash create_local_manifest.sh /vast/projects/bahlo_longstr/nf-str-output/alignments/ input_manifest.tsv 

INPUT_DIR="${1:?Usage: $0 <input_dir> <output_manifest.tsv> [delimiter]}"
OUTPUT_TSV="${2:?Usage: $0 <input_dir> <output_manifest.tsv> [delimiter]}"
DELIMITER="${3:-_}"   # default delimiter is underscore

> "$OUTPUT_TSV"  # truncate/create output file

while IFS= read -r -d '' filepath; do
    basename=$(basename "$filepath")
    name_noext="${basename%.cram}"

    # Extract sample_id: everything before the first delimiter
    sample_id=$(echo "$name_noext" | cut -d"${DELIMITER}" -f1)

    # Determine sequencing type
    type="unknown"
    lower=$(echo "$name_noext" | tr '[:upper:]' '[:lower:]')
    if echo "$lower" | grep -q "ont"; then
        type="ont"
    elif echo "$lower" | grep -qE "pacbio|pb"; then
        type="pacbio"
    elif echo "$lower" | grep -qE "illumina|ilmn"; then
        type="illumina"
    fi

    realpath_file=$(realpath "$filepath")

    printf '%s\t%s\t%s\n' "$sample_id" "$type" "$realpath_file" >> "$OUTPUT_TSV"

done < <(find "$INPUT_DIR" -type f -name "*.cram" -print0 | sort -z)

echo "Manifest written to: $OUTPUT_TSV"

