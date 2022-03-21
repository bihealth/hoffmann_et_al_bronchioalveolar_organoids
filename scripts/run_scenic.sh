#!/bin/bash -e

#arboreto_with_multiprocessing.py HT280.loom hs_hgnc_curated_tfs.txt --method grnboost2 --output HT280_adj.tsv --num_workers 16 --seed 1
pyscenic ctx HT280_adj.tsv *.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl --output HT280_reg.csv --mask_dropouts --num_workers 8 --expression_mtx_fname HT280.loom 
pyscenic aucell HT280.loom HT280_reg.csv --output HT280_pyscenic.loom --num_workers 8
