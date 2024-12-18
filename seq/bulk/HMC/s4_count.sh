#!/bin/bash
featureCounts -p --countReadPairs -t exon -g gene_id -T 10 \
-a /home/hzg/rna/human_genome/gencode.v43.annotation.filtered.egfp.gtf \
-o ./FeatureCounts1.txt BAM/*.sortedByCoord.out.bam
