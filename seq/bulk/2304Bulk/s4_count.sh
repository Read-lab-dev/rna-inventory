#!/bin/bash
featureCounts -p -t exon -g gene_id -T 10 \
-a /home/hzg/backup/lab/human_genome/gencode.v43.annotation.filtered.egfp.gtf \
-o ./FeatureCounts.count BAM/*.sortedByCoord.out.bam
