#!/bin/bash

output_folder="/localdisk/home/s2530615/ICA1/output_group"
count_dir="/localdisk/home/s2530615/ICA1/output_bed/"

process_group() {
    group_file=$1
    output_file="${group_file}_avg.txt"
    
    # 使用 cat 获取所有相关的计数文件内容
    cat $(awk '{gsub(/Tco/,"Tco-"); print "'$count_dir'" $1 "_count.txt"}' $group_file) | 
    # 使用 awk 计算每个基因的总计数、出现次数及保存基因描述
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

