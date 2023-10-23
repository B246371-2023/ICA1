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
