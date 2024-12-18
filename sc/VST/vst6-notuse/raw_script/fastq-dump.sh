cat sample2 | while read id;
do
parallel-fastq-dump --sra-id ${id} --threads 10 --outdir ./ --split-files --gzip --tmpdir /home/hzg/backup/;
done

