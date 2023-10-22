#!/bin/bash
#use bedtoolscoverage to calculate of the reads aligned to the regions specified in the "bedfile"

input_folder="/localdisk/home/s2530615/ICA1/output_bam"
output_folder="/localdisk/home/s2530615/ICA1/output_bed"

mkdir -p "$output_folder"
for file in "$input_folder"/*.bam; do
    if [ -f "$file" ]; then
        filename="${file##*/}"
        filename="${filename%.bam}"
        output_file="$output_folder/${filename}.bed"
        #-a alignment.bam: Specifies the BAM file containing the mapped reads.
        #-b /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed: Specifies the BED file with gene locations.
        #-wa -c: Specifies that you want to write the original entry in the -a file for each overlap and count the number of overlaps. This will give you counts for each read aligned to a gene region.
        bedtools bamtobed -i "$file"  > "$output_file"
        echo "Processed:$filename"
    fi
done
echo "All files processed"
