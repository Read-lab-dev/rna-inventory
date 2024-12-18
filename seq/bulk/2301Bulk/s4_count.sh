#!/bin/bash
featureCounts -p -t exon -g gene_id -T 10 \
-a /mnt/k/Reference/gencode.v43.annotation.gtf \
-o ./FeatureCounts.count BAM/*.sortedByCoord.out.bam
