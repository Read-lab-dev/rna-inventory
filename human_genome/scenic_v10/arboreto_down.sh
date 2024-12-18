#!/bin/bash
tfs=./resources.aertslab.org_cistarget_tf_lists_allTFs_hg38.txt
feather=./hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=./motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
input_loom=./sample.loom
ls $tfs  $feather  $tbl 

#pyscenic step2 cistarget
echo 'Now Processing Cistarget'
pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg_10k.csv \
--num_workers 6  

#pyscenic step3 AUCell
echo 'Now Processing AUCell'
pyscenic aucell \
$input_loom \
reg_10k.csv \
--output out_SCENIC_10k.loom \
--num_workers 10
