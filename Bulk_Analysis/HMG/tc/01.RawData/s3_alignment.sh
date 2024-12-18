#!/bin/bash
mkdir ./BAM
ls *_1.fq.gz | while read id;
do
	echo "*********Processing $id**************"
	STAR --runThreadN 10 --genomeDir /home/hzg/rna/human_genome/index \
         --readFilesCommand zcat --readFilesIn $id ${id%_*}_2.fq.gz \
         --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 \
         --quantMode TranscriptomeSAM GeneCounts \
         --outTmpDir ~/tmp --outFileNamePrefix ../BAM/${id%_*}. \
         --runDirPerm All_RWX
done