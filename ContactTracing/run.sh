WORKDIR=/home/hzg/rna/ContactTracing/   # set this to directory where you keep your input files
mkdir -p $WORKDIR        # create directory if doesn't exist
cd $WORKDIR
git clone https://github.com/LaughneyLab/ContactTracing_tutorial.git
docker pull docker.io/biohpc/scrna2023
docker run --rm -d -v $WORKDIR:/data -p 8888:8888 docker.io/biohpc/scrna2023 /root/scripts/startJupyter.sh /data