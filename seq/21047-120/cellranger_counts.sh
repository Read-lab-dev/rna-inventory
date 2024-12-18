#!/bin/bash
export PATH=/home/hzg/miniconda3/envs/cellranger/cellranger-7.1.0:$PATH
cat sample.txt | while read id;
do
echo "processing $id"
cellranger count --id=$id --fastqs=./ \
  --sample=$id \
  --transcriptome=/home/hzg/rna/human_genome/smallseq/small
done