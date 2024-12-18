#!/bin/bash
mkdir ./BAM
ls *R1.fq.gz | while read id;
do
	echo "*********Processing $id**************"
	STAR --runThreadN 10 --genomeDir /mnt/k/Reference/index \
         --readFilesCommand zcat --readFilesIn $id ${id%_*}_R2.fq.gz \
         --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 \
         --quantMode TranscriptomeSAM GeneCounts \
         --outTmpDir ~/tmp --outFileNamePrefix ./BAM/${id%_*}. \
         --runDirPerm All_RWX
done