#!/bin/bash

# Transfer sam file to bam file
input_folder="/localdisk/home/s2530615/ICA1/output"
output_folder="/localdisk/home/s2530615/ICA1/output_bam"

mkdir -p "$output_folder"

for file in "$input_folder"/*.sam; do
    if [ -f "$file" ]; then
        filename="${file##*/}"
        filename="${filename%.sam}"
        output_file="$output_folder/${filename}.bam"
        sorted_output_file="$output_folder/${filename}_sorted.bam"

        # Convert SAM to BAM
        samtools view -b "$file" > "$output_file"

        # Sort BAM
        samtools sort "$output_file" -o "$output_file"

        # Index sorted BAM
        samtools index "$output_file"

        echo "Processed: $filename"
    fi
done

echo "All files processed"
