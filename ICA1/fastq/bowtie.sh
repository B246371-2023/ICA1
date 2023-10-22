for file in $/localdisk/home/s2530615/ICA1/fastq;do
    if [ -f "$file" ]; then
        filename="${file##*.}"
        filename="${filename%_1.fq}"
        bowtie -x Tcongo_IL3000 -1 "${folder_path}/${filename}_1.fq" -2 "${folder_path}/${filename}_2.fq" -S "${folder_path}/${filename}.bam"
        echo "Processed: $filename"
    fi
done
