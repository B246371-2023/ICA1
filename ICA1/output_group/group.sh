# Input file
input_file="/localdisk/home/s2530615/ICA1/fastq/Tco2.fqfiles"
output_folder="/localdisk/home/s2530615/ICA1/output_group"
mkdir -p "$output_folder"

# Function to extract data based on criteria
filter_data() {
    sample_type=$1
    time=$2
    treatment=$3
    output_file="$output_folder/${sample_type}_${time}h_${treatment}.txt"

    # Extract rows matching criteria
    awk -v sample="$sample_type" -v time="$time" -v treatment="$treatment" \
        'BEGIN{OFS="\t"} $2 == sample && $4 == time && $5 == treatment' \
        "$input_file" > "$output_file"

    # If the resulting file is empty, remove it
    if [[ ! -s "$output_file" ]]; then
        rm -rf "$output_file"
        echo "No data for ${sample_type}_${time}h_${treatment}"
    else
        echo "Data saved to $output_file"
    fi
}

# Filtering data for each group combination and saving to separate files
for sample in "WT" "Clone1" "Clone2"; do
    for time in "0" "24" "48"; do
        for treatment in "Uninduced" "Induced"; do
            filter_data "$sample" "$time" "$treatment"
        done
    done
done

