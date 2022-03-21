# Code repository for Hoffmann et al. "Human alveolar progenitors generate dual lineage bronchioalveolar organoids"

## data access

processed data is available from NCBI GEO under accession [GSE197949](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197949)  

- scRNAseq: cellranger output (raw h5 files) for lung organoids (GSM5934385..406)
- bulk RNAseq: counts for lung organoids (GSM5934407..14)

and should be downloaded and stored in the `cellranger` and `feature_counts` subfolders like so (GSM[0-9]*_ accession codes should be stripped off!)

```
data
│── cellranger
│   ├── AO13_pool_batch1_control_raw.h5
│   ├── AO13_pool_batch2_control_raw.h5
│   └── ...
└── feature_counts
    ├── ctrl_CHIR_1.feature_counts
    ├── ctrl_CHIR_2.feature_counts
    └── ...
```

raw h5 files (`lung_{102C,219V}_{BSA,H3N2}_raw.h5`, GSM59582{61,63,77,79}) for lung tissue explants are available [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198864)

final R objects are also available and should go into `data/seurat`

## CellBender

[CellBender](https://github.com/broadinstitute/CellBender) `remove-background` is used to remove viral RNA background especially from autopsy single-nuc samples

it could be run like so:

```
while read -r sample ncells; do
    bash scripts/run_cellbender.sh ${sample} ${ncells} 
done < data/cellranger/expected_cells.txt
```

## Seurat processing

R code in `R` is split into several markdowns that can be simply knitted:

- `process_HT280_Epcam.Rmd`
- `process_organoid_all.Rmd`
- `process_organoid_control.Rmd`
- `process_organoid_HT2.Rmd`
- `process_organoid_pool.Rmd`
- `process_lung_tissue.Rmd`

## SCENIC

[SCENIC](https://github.com/aertslab/SCENIC) is run on HT280+/Epcam+ sorted cells

required databases:

- hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
- hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
- hs_hgnc_curated_tfs.txt 
- motifs-v9-nr.hgnc-m0.001-o0.0.tbl

1. prepare loom input files `python scripts/prepare_scenic.py`
2. run scenic: `cd data/scenic; bash ../../scripts/run_scenic.sh`
3. process results: `python scripts/process_scenic.py`

## CONICS

[CONICS](https://github.com/diazlab/CONICS) is used to check for potential copy-number variation events

required input files:

- chromosome_arm_positions_grch38.txt
- gene_positions.txt  

## external datasets

lung cell atlas reference

- `droplet_normal_lung_seurat_ntiss10x.P2.anno.20191002.RC4.Robj` from [here](https://www.synapse.org/#!Synapse:syn21560412) converted to Seurat V3 and stored as RDS

## paper figures

all code for paper figures is in `paper_figures.Rmd`
