cat sample | while read id;
do
cellranger count --id=$id --fastqs=./ \
  --sample=$id \
  --transcriptome=/home/hzg/backup/lab/human_genome/homo_sapiens_h3
done