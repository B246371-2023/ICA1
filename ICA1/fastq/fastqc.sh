#fastqc.sh
total_files=0
all_passed_files=0
for zipfile in Tco-*_fastqc.zip; do
   
    tmpdir=$(mktemp -d)

    unzip -q "$zipfile" -d "$tmpdir"

    datafile="$tmpdir/$(basename "$zipfile" .zip)"/fastqc_data.txt
    summaryfile="$tmpdir/$(basename "$zipfile" .zip)"/summary.txt

    outputfile="$(basename "$zipfile" _fastqc.zip)_access.txt"

    head -10 "$datafile" > "$outputfile"

    if [[ -f "$summaryfile" ]]; then
        cat "$summaryfile" >> "$outputfile"
    fi
    if ! grep -q "^FAIL" "$summaryfile" && ! grep -q "^WARN" "$summaryfile"; then
        ((all_passed_files++))
    fi
    ((total_files++))

    rm -rf "$tmpdir"
done


echo "Total number of files: $total_files"
if [ $all_passed_files -eq 0 ]; then
    echo "Every sequence quality is good for all files."
fi
rm -rf *.zip
