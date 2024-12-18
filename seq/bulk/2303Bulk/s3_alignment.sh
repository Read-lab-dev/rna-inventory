#!/bin/bash
cat cat.txt | while read id;
do
	echo "*********Processing $id**************"
	STAR --runThreadN 10 --genomeDir /mnt/k/Reference/index \
         --readFilesCommand zcat --readFilesIn $id ${id%_*}_R2.fq.gz \
         --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 \
         --quantMode TranscriptomeSAM GeneCounts \
         --outTmpDir ~/tmp --outFileNamePrefix ${id%_*}. \
         --runDirPerm All_RWX
done

featureCounts -p -t exon -g gene_id -T 10 \
-a /mnt/k/Reference/gencode.v43.annotation.gtf \
-o ./FeatureCounts.count fastp/*.sortedByCoord.out.bam