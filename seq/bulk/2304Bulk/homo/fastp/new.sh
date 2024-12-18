ls *R1.fq.gz | while read id;
do
	echo "*********Processing $id**************"
	STAR --runThreadN 10 --genomeDir /home/hzg/backup/lab/human_genome/STAR \
         --readFilesCommand zcat --readFilesIn $id ${id%_*}_R2.fq.gz \
         --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 \
         --quantMode TranscriptomeSAM GeneCounts \
         --outTmpDir /home/hzg/backup/tmp --outFileNamePrefix ./BAM/${id%_*}. \
         --runDirPerm All_RWX
done

featureCounts -p -t exon -g gene_id -T 10 \
-a /home/hzg/backup/lab/human_genome/gencode.v43.annotation.filtered.egfp.gtf \
-o ./FeatureCounts.count BAM/*.sortedByCoord.out.bam
