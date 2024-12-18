#!/bin/bash
featureCounts -p -t exon -g gene_id -T 10 \
-a /mnt/k/Reference/gencode.v43.annotation.gtf \
-o /mnt/c/Users/read_lab/"OneDrive - Emory University"/bulk/integrated/FeatureCounts.count *.sortedByCoord.out.bam
