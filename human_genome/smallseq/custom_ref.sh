cellranger mkgtf \
  gencode.v43.annotation.gtf \
  gencode.v43.annotation.filtered.gtf \
  --attribute=gene_biotype:protein_coding

 cellranger mkref \
  --genome=homo_sapiens \
  --fasta=GRCh38.primary_assembly.genome.fa \
  --genes=gencode.v43.annotation.filtered.gtf