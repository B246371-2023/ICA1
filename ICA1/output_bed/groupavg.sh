#!/bin/bash

samples_info_path="/localdisk/home/s2530615/ICA1/fastq/Tco2.fqfiles"
count_dir="/localdisk/home/s2530615/ICA1/output_bed/"

while IFS=$'\t' read -r SampleName SampleType Replicate Time Treatment End1 End2; do
    
    if [ "$SampleName" != "SampleName" ]; then
        group="${SampleType}_${Replicate}_${Time}_${Treatment}"
        output_file="${group}_average.txt"
        
     
        > "$output_file"
        
       
        modified_sample_name="Tco-${SampleName#Tco}"
        count_file="${count_dir}${modified_sample_name}_counts.txt"
        
        if [ -f "$count_file" ]; then
           
            awk -v group="$group" '{
                gene_name=$4;
                gene_description=$5;
                count=$6;
                gene_totals[gene_name] += count;
                gene_counts[gene_name]++;
                gene_desc[gene_name] = gene_description;
            } END {
                for (gene in gene_totals) {
                    avg = gene_totals[gene] / gene_counts[gene];
                    printf "%s\t%s\t%.2f\n", gene, gene_desc[gene], avg >> group"_average.txt";
                }
            }' "$count_file"
        else
            echo "Count file $count_file not found!"
        fi
    fi
done < "$samples_info_path"

