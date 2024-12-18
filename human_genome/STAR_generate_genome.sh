#!/bin/bash
fasta=/home/hzg/backup/lab/human_genome/GRCh38.primary_assembly.genome.egfp.fa
gtf=/home/hzg/backup/lab/human_genome/gencode.v43.annotation.filtered.egfp.gtf

# Generate index Of Brassica napus

STAR --runMode genomeGenerate \
--runThreadN 10 \
--genomeFastaFiles $fasta \
--genomeDir ./STAR \
--sjdbGTFfile $gtf \
--genomeSAindexNbases 14 \
--sjdbOverhang 99

