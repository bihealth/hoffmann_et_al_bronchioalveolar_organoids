import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import loompy as lp

cells=[line.split()[0] for line in open('data/seurat/HT280_Epcam_barcodes.txt')]
anno=pd.read_csv('data/seurat/HT280_Epcam_meta.csv',header=0,index_col=0)
samples=anno['orig.ident'].unique()

samples=dict(zip(['AO22_HT280','AO23_HT280','AO34_HT280_Epcam',
                  'AO35_HT280_Epcam','AO36_HT280_Epcam','AO37_HT280_Epcam'],samples))

ads=[]
for lib,lib2 in samples.items():
    h5_file=os.path.join('data','cellbender',lib+'_filtered.h5')
    if not os.path.isfile(h5_file):
        print(h5_file)
        continue
    ad=sc.read_10x_h5(h5_file)
    ad.X=scipy.sparse.csr_matrix(ad.X)
    ad.obs.index=lib2.replace('_','-').replace('HT280-Epcam','HT280Epcam')+'_'+ad.obs.index.astype(str)
    ad.obs['sample']=lib2
    ad.var_names_make_unique()
    ads.append(ad)
ad=ads[0].concatenate(ads[1:],batch_key=None,index_unique=None)

ad=ad[cells]
row_attrs = {
    "Gene": np.array(ad.var_names) ,
}
col_attrs = {
    "CellID": np.array(ad.obs_names) ,
    "sample": np.array(ad.obs['sample']),
    "nGene": np.array( np.sum(ad.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(ad.X.transpose() , axis=0)).flatten() ,
}
lp.create('data/scenic/HT280.loom', ad.X.transpose(), row_attrs, col_attrs)

