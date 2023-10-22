#!/bin/bash

# 前提设置
output_folder="/localdisk/home/s2530615/ICA1/output_group"
count_dir="/localdisk/home/s2530615/ICA1/output_bed/"

process_group() {
    group_file=$1

    # 在此处，每次函数开始时重置这两个数组
    declare -A gene_totals
    declare -A gene_counts

    echo "==== Processing group file: $group_file ===="

    # 检查文件是否存在并可读
    if [[ ! -r "$group_file" ]]; then
        echo "Cannot read group file: $group_file"
        return
    fi

    while IFS=$'\t' read -r -a fields; do
        sample="${fields[0]}"
        modified_sample_name="${sample/Tco/Tco-}"
        count_file="${count_dir}${modified_sample_name}_count.txt"
        
        echo "  -- Processing count file for sample: $sample --"

        if [[ -f "$count_file" ]]; then
            # 使用管道代替临时文件
            cut -f 4,6 "$count_file" | awk 'BEGIN{OFS="\t"} {print $1, $2}' | while IFS=$'\t' read -r gene count; do
                gene_totals["$gene"]=$(( gene_totals["$gene"] + count ))
                gene_counts["$gene"]=$(( gene_counts["$gene"] + 1 ))
            done
        else
            echo "Warning: Count file does not exist for sample: $sample"
        fi
    done < "$group_file"

    # 将当前组的总结写入一个特定的文件中
    output_file="${group_file%.txt}_avg.txt"
    echo "  -- Writing output to: $output_file --"

    for gene in $(echo "${!gene_totals[@]}" | tr ' ' '\n' | sort); do
        total="${gene_totals[$gene]}"
        count="${gene_counts[$gene]}"
        average=$(echo "scale=2; $total / $count" | bc)
        echo -e "$gene\t$count\t$total\t$average" >> "$output_file"
    done
    echo "==== Finished processing group file: $group_file ===="
    echo
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

