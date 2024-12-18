cat sample | while read id;
do
velocyto run10x -m ~/rna/human_genome/hg38_mask.gtf /home/hzg/rna/sequencing/21047FL-89/$id ~/rna/human_genome/gencode.v43.annotation.filtered.egfp.gtf -@ 10 -v -t uint32
done