#!/bin/bash
mkdir ./BAM
STAR --genomeLoad LoadAndExit --genomeDir /home/hzg/rna/human_genome/STAR
ls *R1.fq.gz | while read id;
do
	echo "*********Processing $id**************"
	STAR --runThreadN 10 --genomeDir /home/hzg/rna/human_genome/STAR \
         --readFilesCommand zcat --readFilesIn $id ${id%_*}_R2.fq.gz \
         --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 \
         --quantMode TranscriptomeSAM GeneCounts \
         --outTmpDir /home/hzg/backup/tmp --outFileNamePrefix ./BAM/${id%_*}. \
         --runDirPerm All_RWX
done
STAR --genomeLoad Remove --genomeDir /home/hzg/rna/human_genome/STAR
