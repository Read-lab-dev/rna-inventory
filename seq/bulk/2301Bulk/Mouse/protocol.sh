#!/bin/bash
ls *_R1.fastq.gz | while read id;
 do 
    echo "***********Now processing $id****************"
    fastp -i $id -o ./${id%_*}_R1.fq.gz -I ${id%_*}_R2.fastq.gz -O ./${id%_*}_R2.fq.gz -h ${id%_*}.html -w 10;
 done

rm -r ~/tmp
ls *_R1.fq.gz | while read id;
do
   echo "*********Processing $id**************"
   STAR --runThreadN 10 --genomeDir ~/rna/mouse_genome/index \
        --readFilesCommand zcat --readFilesIn $id ${id%_*}_R2.fq.gz \
        --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 \
        --outTmpDir ~/tmp --outFileNamePrefix ./${id%_*}. \
        --runDirPerm All_RWX
done

featureCounts -p -t exon -g gene_id -T 10 \
-a ~/rna/mouse_genome/gencode.vM33.annotation.gtf \
-o ./FeatureCounts.count *.sortedByCoord.out.bam
