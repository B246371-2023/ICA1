#!/bin/bash

# Create necessary directories for the workflow
mkdir -p ./ICA1
mkdir -p ./ICA1/fastq

# Copy fastq data from the specified source. If copying fails, exit the script with an error.
cp -r /localdisk/data/BPSM/ICA1/fastq ./ICA1 || { echo "Failed to copy files"; exit 1; }

# Notify the user about the start of decompression
echo "Decompressing files..."

# Decompress the fastq files and check for errors
gunzip ./ICA1/fastq/*.fq.gz  || { echo "Failed to decompress fastq files"; exit 1; }

# Notify the user that decompression is complete
echo "Decompression completed"

# Analyze the quality of the decompressed sequences using FastQC
fastqc ./ICA1/fastq/*.fq


# fastqc.sh section starts here
# Initialize counters for total files and files that passed FastQC
total_files=0
all_passed_files=0

# Loop through each FastQC output zip file
for zipfile in Tco-*_fastqc.zip; do
    # Create a temporary directory for extracting contents of the zip file
    tmpdir=$(mktemp -d)
    
    # Extract the contents of the zip file to the temporary directory
    unzip -q "$zipfile" -d "$tmpdir"
    
    # Define the path to the fastqc data and summary files within the extracted folder
    datafile="$tmpdir/$(basename "$zipfile" .zip)"/fastqc_data.txt
    summaryfile="$tmpdir/$(basename "$zipfile" .zip)"/summary.txt
    outputfile="$(basename "$zipfile" _fastqc.zip)_access.txt"
    
    # Notify the user about the current file being processed
    echo "Processing: $zipfile"  
    
    # Copy the first 10 lines of the data file to the output file
    head -10 "$datafile" > "$outputfile"
    
    # If the summary file exists, append its contents to the output file
    if [[ -f "$summaryfile" ]]; then
        cat "$summaryfile" >> "$outputfile"
    fi
    
    # Check the summary file for any "FAIL" or "WARN" status, if none found increment the passed files counter
    if ! grep -q "^FAIL" "$summaryfile" && ! grep -q "^WARN" "$summaryfile"; then
        ((all_passed_files++))
    fi
    
    # Increment the total files counter
    ((total_files++))
    
    # Delete the temporary directory
    rm -rf "$tmpdir"
done

# Display the total number of processed files
echo "Total number of files: $total_files"

# If all files passed the FastQC check, notify the user
if [ $all_passed_files -eq 0 ]; then
    echo "Every sequence quality is good for all files."
fi

# Delete all the zip files
rm -rf *.zip


# Copy the genome file to the target directory
cp /localdisk/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz ./ICA1/fastq/Tcongo_IL3000.fasta.gz  || { echo "Failed to copy genome files"; exit 1; }
# Decompress the genome file
gunzip ./ICA1/fastq/Tcongo_IL3000.fasta.gz

# bowtie2.sh section
# Index the genome using bowtie2 for alignment purposes
bowtie2-build ./ICA1/fastq/Tcongo_IL3000.fasta ./ICA1/fastq/Tcongo_IL3000

# Define variables for paths and computational resources
index="./ICA1/fastq/Tcongo_IL3000"
input_folder="./ICA1/fastq"
output_folder="./ICA1/output"
# Get the number of threads (cores) to be used
threads=256

# Create output directory if not existing
mkdir -p "$output_folder"

# Loop through each pair of fastq files for alignment
for file in "$input_folder"/*_1.fq; do
    if [ -f "$file" ]; then
        filename="${file##*/}"
        filename="${filename%_1.fq}"
        
        output_file="$output_folder/${filename}.sam"
        
        # Align the paired-end reads to the indexed genome using bowtie2
        bowtie2 -x "$index" -1 "$file" -2 "${file%_1.fq}_2.fq" -S "$output_file" -p $threads || { echo "Failed to align with bowtie2 for $filename"; exit 1; }
        
        echo "Processed: $filename"
    fi
done

echo "All files processed bowtie2 to align."

# samtools2 section
# Convert and process alignment files
input_folder="./ICA1/output"
output_folder="./ICA1/output_bam"
threads=128
mkdir -p "$output_folder"

# Loop through each SAM file for conversion and processing
for file in "$input_folder"/*.sam; do
    if [ -f "$file" ]; then
        filename="${file##*/}"
        filename="${filename%.sam}"
        output_file="$output_folder/${filename}.bam"
        sorted_output_file="$output_folder/${filename}.bam"

        # Convert SAM to BAM
        samtools view -b -@ $threads "$file" > "$output_file"
        
        # Sort BAM files
        samtools sort -@ $threads "$output_file" -o "$sorted_output_file"

        # Index the sorted BAM files
        samtools index "$sorted_output_file" 
        
        echo "Processed: $filename"
    fi
done

echo "All files processed samtools to bam."

# Copy BED annotation file to the target directory
cp /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed ./ICA1 || { echo "Failed to copy bed file"; exit 1; }


# bedtools4.sh
# Use bedtools to count reads aligned to specific regions specified in a BED file.

# Define the bed file and folders for input and output
bedfile="./ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed"
input_folder="./ICA1/output_bam"
output_folder="./ICA1/output_bed"

# Create the output directory
mkdir -p "$output_folder"

# Loop through BAM files and use bedtools to count reads overlapping regions in the bed file
for file in "$input_folder"/*.bam; do
    if [ -f "$file" ]; then
        filename="${file##*/}"
        filename="${filename%.bam}"
        output_file="$output_folder/${filename}_counts.txt"
    
        bedtools multicov -bams "$file" -bed "$bedfile" > "$output_file"
        echo "Processed:$filename"
    fi
done

echo "All files processed with bedtools for counting"


#filter group
# Group data based on sample type, time, and treatment
input_file="./ICA1/fastq/Tco2.fqfiles"
output_folder="./ICA1/output_group"
mkdir -p "$output_folder"

# Function to filter and extract data based on sample type, time, and treatment
filter_data() {
    # Extract parameters for filtering
    sample_type=$1
    time=$2
    treatment=$3
    output_file="$output_folder/${sample_type}_${time}h_${treatment}.txt"
    
    # Use awk to filter and extract the rows matching criteria
    awk -v sample="$sample_type" -v time="$time" -v treatment="$treatment" \
    'BEGIN{FS=OFS="\t"} $2 == sample && $4 == time && $5 == treatment' \
    $input_file > $output_file 
    
    # If the resulting file is empty, remove it and log the message
    if [[ ! -s $output_file ]]; then
        rm -rf "$output_file"
        echo "No data for ${sample_type}_${time}h_${treatment}"
    else
        echo "Data saved to $output_file"
    fi
}

# Loop through each sample type, time, and treatment combination to filter data
for sample in "WT" "Clone1" "Clone2"; do
    for time in "0" "24" "48"; do
        for treatment in "Uninduced" "Induced"; do
            filter_data "$sample" "$time" "$treatment"
        done
    done
done

# group mean
# Calculate mean of counts for each gene in grouped data
output_folder="./ICA1/output_group"
count_dir="./ICA1/output_bed/"
mkdir -p "$output_folder"

# Function to process each group and calculate the average count for each gene
process_group() {
    group_file=$1
    output_file="${group_file%.txt}_avg.txt"
    
    # Concatenate count files, calculate and save average counts for each gene
    cat $(awk '{gsub(/Tco/,"Tco-"); print "'$count_dir'" $1 "_counts.txt"}' $group_file) | 
    # Include cut gene describtion by loop 
    awk '{
        gene = $4;
        description = "";
        for(i=5; i<NF; i++) {
            description = description $i " ";
        }
        count = $NF;
        total[gene] += count;
        count_gene[gene]++;
        gene_desc[gene] = description;
    } END {
        for (gene in total) {
            average = total[gene] / count_gene[gene];
            print gene, gene_desc[gene], average;
        }
    }' > "$output_file"
}

# Loop through each sample type, time, and treatment combination to process grouped data
for sample_type in "WT" "Clone1" "Clone2"; do
    for time in "0" "24" "48"; do
        for treatment in "Uninduced" "Induced"; do
            group_file="${output_folder}/${sample_type}_${time}h_${treatment}.txt"
            if [[ -f "$group_file" ]]; then
                process_group "$group_file"
            else
                echo "Warning: Group file does not exist: $group_file"
            fi
        done
    done
done


# fold-change
# Calculate the fold changes between all possible pairs of conditions and store the results.

# Define the location where the grouped data is stored
output_folder="/localdisk/home/s2530615/ICA1/output_group"

# Nested loops to go through all possible combinations of conditions to compare 
for sample_type in "WT" "Clone1" "Clone2"; do
    for time in "0" "24" "48"; do
        for treatment in "Uninduced" "Induced"; do
            # Define the first file based on the current iteration's conditions
            file1="${output_folder}/${sample_type}_${time}h_${treatment}_avg.txt"
            
            for sample_type2 in "WT" "Clone1" "Clone2"; do
                for time2 in "0" "24" "48"; do
                    for treatment2 in "Uninduced" "Induced"; do
                        # Define the second file based on the current iteration of the inner loops
                        file2="${output_folder}/${sample_type2}_${time2}h_${treatment2}_avg.txt"
                        
                        # Make sure both files exist and aren't the same
                        if [[ -f "$file1" && -f "$file2" && "$file1" != "$file2" ]]; then
                            # Define the output file name based on the two files being compared
                            output_file="${output_folder}/fold_changes_${sample_type}_${time}h_${treatment}_vs_${sample_type2}_${time2}h_${treatment2}.txt"
                            echo "Processing: $file1 vs $file2"
                            
                            # Use awk to compute fold changes. Read each line of file1, then get the matching line from file2.
                            # Calculate fold change ratio and include gene description.
                            # If denominator is zero, set fold change ratio to "inf" (infinity)
                            awk '{
                                gene = $1;
                                description = "";
                                for (i = 2; i < NF; i++) {
                                    description = description $i " ";
                                }
                                fold_change = $NF;
                                getline < "'$file2'";
                                fold_change2 = $NF;
                                fold_change_ratio = (fold_change2 != 0) ? fold_change / fold_change2 : "inf";
                                print gene, fold_change_ratio, description;
                            }' "$file1" | sort -k2,2nr > "$output_file"  # Sorting the output by fold change in descending order

                            echo "Processing complete: $file1 vs $file2"
                        fi
                    done
                done
            done
        done
    done
done

