//bedtools.sh
#!/bin/bash
#use bedtoolscoverage to calculate of the reads aligned to the regions specified in the "bedfile"
bedfile="/localdisk/home/s2530615/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed"
input_folder="/localdisk/home/s2530615/ICA1/output_bam"
output_folder="/localdisk/home/s2530615/ICA1/output_bed"

mkdir -p "$output_folder"
for file in "$input_folder"/*.bam; do
    if [ -f "$file" ]; then
        filename="${file##*/}"
        filename="${filename%.bam}"
        output_file="$output_folder/${filename}_count.txt"
    
        bedtools intersect -a "$bedfile" -b "$file" -c > "$output_file"
        echo "Processed:$filename"
    fi
done
echo "All files processed"
