touch H3K27M_allid.txt
touch H3K27WT_allid.txt
cat sample | while read id;
do
  samtools view  $id/outs/possorted_genome_bam.bam chr1:226064352-226064479 | grep CGCATGAGT |grep CB:Z: |grep -Eo 'CB:Z:[A-Z,1-9]{1,18}'| sed 's/\CB:Z://g'|sed 's/$/-1/' >>H3K27M_allid.txt
  samtools view  $id/outs/possorted_genome_bam.bam chr1:226064352-226064479 | grep CGCAAGAGT |grep CB:Z: |grep -Eo 'CB:Z:[A-Z,1-9]{1,18}'| sed 's/\CB:Z://g'|sed 's/$/-1/' >>H3K27WT_allid.txt
done

 #rm ${id}_mut3.3.id.txt
 #rm ${id}_wild3.3.id.txt