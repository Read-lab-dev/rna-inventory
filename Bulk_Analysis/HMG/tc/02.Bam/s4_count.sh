#!/bin/bash
featureCounts -p -t exon -g gene_id -T 10 \
-a /home/hzg/rna/human_genome/gencode.v43.annotation.gtf \
-o ./FeatureCounts.count *.bam
