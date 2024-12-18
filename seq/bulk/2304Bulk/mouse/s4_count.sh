#!/bin/bash
featureCounts -p -t exon -g gene_id -T 10 \
-a /mnt/k/mouse_genome/gencode.vM33.annotation.gtf \
-o ./FeatureCounts.count BAM/*.sortedByCoord.out.bam
