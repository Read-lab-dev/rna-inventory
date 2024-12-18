cat sample | while read id ;
do
 echo https://sra-pub-run-odp.s3.amazonaws.com/sra/${id}/${id} >>link
done