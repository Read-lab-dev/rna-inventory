#!/bin/bash
fasta=/mnt/e/Bulk/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=/mnt/e/Bulk/Reference/Homo_sapiens.GRCh38.99.gtf

# Generate index Of Brassica napus

STAR --runMode genomeGenerate \
--runThreadN 40 \
--genomeFastaFiles $fasta \
--genomeDir /mnt/e/Bulk/Reference/index \
--sjdbGTFfile $gtf \
--genomeSAindexNbases 14 \
--sjdbOverhang 99

