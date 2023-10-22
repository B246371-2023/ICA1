#Index the genome for bowtie2 
bowtie2-build /home/s2530615/ICA1/fastq/Tcongo_IL3000.fasta Tcongo_IL3000

index="Tcongo_IL3000"

input_folder="/localdisk/home/s2530615/ICA1/fastq"

output_folder="/localdisk/home/s2530615/ICA1/output"

#lscpu look the whole system cpu information and find the core cpu number and
threads=256

mkdir -p "$output_folder"

for file in "$input_folder"/*_1.fq; do
    if [ -f "$file" ]; then
        filename="${file##*/}"
        filename="${filename%_1.fq}"
        
        output_file="$output_folder/${filename}.sam"
        
        #bowtie2 align the read pairs to the T.C.
        #we know it a paired-end sequence so we use parameter -x -1 -2
        bowtie2 -x "$index" -1 "$file" -2 "${file%_1.fq}_2.fq" -S "$output_file" -p $threads || { echo "Failed to align with bowtie2 for $filename"; exit 1; }
        
        echo "Processed: $filename"
    fi
done

echo "All files processed bowtie2 to align."
~
~
~
~
~
~
~
~
~
~
~
~

