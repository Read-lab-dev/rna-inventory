 #!/bin/bash

 samples=("CzL" "CzS" "DMSO")
# Define paths and variables
for index in {1..3}; do
    SAMPLE="${samples[index-1]}"
    index_str=$(printf "%02d" $index)
    BAM="/home/hzg/backup/inhouse_SC_seq/21047FL-120/21047FL-120-${index_str}-01-01-$SAMPLE/outs/possorted_genome_bam.bam"
    samtools view  $BAM \
 chr7:116739950-122004660 | grep ATAATCAAGGTTGTGGTT |grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c|awk 'BEGIN {OFS="\t"} {print $0, "'$SAMPLE'"}' > reads_per_barcode_${SAMPLE}_ZM.txt
done
cat *_ZM.txt >reads_per_barcode_ZM.txt

INPUT_FILE="reads_per_barcode_ZM.txt"
OUTPUT_FILE="reads_per_barcode_ZM-fusion.txt"

echo -e "Count\tBarcode\tSample" > "$OUTPUT_FILE"

awk '{print $1 "\t" $2 "\t" $3}' "$INPUT_FILE" >> "$OUTPUT_FILE"

echo "Success: $OUTPUT_FILE"

rm *_ZM.txt
