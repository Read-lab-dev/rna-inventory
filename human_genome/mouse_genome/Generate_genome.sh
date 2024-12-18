#!/bin/bash
fasta=/mnt/k/mouse_genome/GRCm39.primary_assembly.genome.fa
gtf=/mnt/k/mouse_genome/gencode.vM33.annotation.gtf

# Generate index Of Brassica napus

STAR --runMode genomeGenerate \
--runThreadN 10 \
--genomeFastaFiles $fasta \
--genomeDir /mnt/k/mouse_genome/index \
--sjdbGTFfile $gtf \
--genomeSAindexNbases 14 \
--sjdbOverhang 99

