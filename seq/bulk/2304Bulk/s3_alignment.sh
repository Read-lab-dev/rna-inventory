#!/bin/bash
mkdir ./BAM
ls ./fastp/*R1.fq.gz | while read id;
do
	echo "*********Processing $id**************"
	STAR --runThreadN 10 --genomeDir /home/hzg/backup/lab/human_genome/STAR \
         --readFilesCommand zcat --readFilesIn $id ${id%_*}_R2.fq.gz \
         --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 \
         --quantMode TranscriptomeSAM GeneCounts \
         --outTmpDir /home/hzg/backup/tmp --outFileNamePrefix ./BAM/${id%_*}. \
         --runDirPerm All_RWX
done