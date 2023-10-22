#!/bin/bash

# 前提设置
output_folder="/localdisk/home/s2530615/ICA1/output_group"
count_dir="/localdisk/home/s2530615/ICA1/output_bed/"
bed_file="/localdisk/home/s2530615/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed"

process_group() {
    group_file="$1"
    # 每次函数开始时重置这两个数组
    declare -A gene_totals
    declare -A gene_counts

    echo "Processing group file: $group_file"

    # 读取分组文件中的每个样本并处理
    while IFS=$'\t' read -r -a fields; do
        sample="${fields[0]}"
        modified_sample_name="${sample/Tco/Tco-}"
        count_file="${count_dir}${modified_sample_name}_count.txt"
        
        echo "Processing count file for sample: $sample"

        if [[ -f "$count_file" ]]; then
            awk 'BEGIN{OFS="\t"} {
                gene = $1;
                count = $2;
                print gene, count;
            }' <(cut -f 4,6 "$count_file") > temp_counts.txt

            # 更新局部数组
            while IFS=$'\t' read -r gene count; do
                gene_totals["$gene"]=$(( gene_totals["$gene"] + count ))
                gene_counts["$gene"]=$(( gene_counts["$gene"] + 1 ))
            done < temp_counts.txt
            rm -f temp_counts.txt
        else
            echo "Warning: Count file does not exist for sample: $sample"
        fi
    done < "$group_file"

    # 将当前组的总结写入一个特定的文件中，并进行排序
    output_file="${group_file%.txt}_avg.txt"
    for gene in $(echo "${!gene_totals[@]}" | tr ' ' '\n' | sort); do
        total="${gene_totals[$gene]}"
        count="${gene_counts[$gene]}"
        average=$(echo "scale=2; $total / $count" | bc)
        echo -e "$bed_line\t$total\t$count\t$average" >> "$output_file.tmp"
    done
    # 对输出文件排序，然后加上.bed文件的第5列，并覆盖原始输出文件
    sort "$output_file.tmp" | awk -v bed_file="$bed_file" 'BEGIN{OFS="\t"} {
        bed_line=$(cut -f 5 bed_file);
        print $0, bed_line;
    }' > "$output_file"
    rm "$output_file.tmp"
}

# 对每个样本类型、时间和处理情况进行处理
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
