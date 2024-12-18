#!/bin/bash
ls ./homo/*_R1.fastq.gz | while read id;
 do 
    echo "***********Now processing $id****************"
    fastp -i $id -o ./fastp/${id%_*}_R1.fq.gz -I ${id%_*}_R2.fastq.gz -O ./fastp/${id%_*}_R2.fq.gz -h ${id%_*}.html -w 10;
 done