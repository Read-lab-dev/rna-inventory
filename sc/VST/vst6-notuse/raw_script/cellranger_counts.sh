export PATH=/home/hzg/miniconda3/envs/cellranger/cellranger-7.1.0:$PATH
cat sample | while read id;
do
echo "processing $id"
cellranger count --id=$id --fastqs=./ \
  --sample=$id \
  --transcriptome=/home/hzg/rna/human_genome/homo_sapiens_egfp
done