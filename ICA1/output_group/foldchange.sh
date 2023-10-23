#!/bin/bash

# File storage location
output_folder="/localdisk/home/s2530615/ICA1/output_group"

# Loop to compare and calculate fold changes based on different conditions
for sample_type in "WT" "Clone1" "Clone2"; do
    for time in "0" "24" "48"; do
        for treatment in "Uninduced" "Induced"; do
            file1="${output_folder}/${sample_type}_${time}h_${treatment}_avg.txt"
            for sample_type2 in "WT" "Clone1" "Clone2"; do
                for time2 in "0" "24" "48"; do
                    for treatment2 in "Uninduced" "Induced"; do
                        file2="${output_folder}/${sample_type2}_${time2}h_${treatment2}_avg.txt"
                        # Check if both files exist
                        if [[ -f "$file1" && -f "$file2" && "$file1" != "$file2" ]]; then
                            # Calculate fold change and output to a file
                            output_file="${output_folder}/fold_changes_${sample_type}_${time}h_${treatment}_vs_${sample_type2}_${time2}h_${treatment2}.txt"
                            echo "Processing: $file1 vs $file2"
                            # Use awk to calculate fold change and include full gene description
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
                            }' "$file1" | sort -k2,2nr > "$output_file"  # Sort by the second column in descending order
                            echo "Processing complete: $file1 vs $file2"
                        fi
                    done
                done
            done
        done
    done
done

