import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import loompy as lp
import json
import zlib
import base64
from scipy.spatial.distance import jensenshannon
from math import sqrt

anno=pd.read_csv('data/seurat/HT280_Epcam_meta.csv',header=0,index_col=0)

# collect SCENIC AUCell output
lf = lp.connect('data/scenic/HT280_pyscenic.loom', mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
lf.close()

auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
    tmp = x.get('regulon').replace("(","_(")
    x.update( {'regulon': tmp} )

auc_mtx.to_csv('data/scenic/HT280_scenic_AUC.csv')
pd.DataFrame(regulons).to_csv('data/scenic/HT280_scenic_regulons.csv')

def regulon_specificity_scores(auc_mtx, cell_type_series):
    """
    Calculates the Regulon Specificty Scores (RSS). [doi: 10.1016/j.celrep.2018.10.045]
    :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
    :param cell_type_series: A pandas Series object with cell identifiers as index and cell type labels as values.
    :return: A pandas dataframe with the RSS values (cell type x regulon).
    """

    cell_types = list(cell_type_series.unique())
    n_types = len(cell_types)
    regulons = list(auc_mtx.columns)
    n_regulons = len(regulons)
    rss_values = np.empty(shape=(n_types, n_regulons), dtype=np.float)

    def rss(aucs, labels):
        # jensenshannon function provides distance which is the sqrt of the JS divergence.
        return 1.0 - jensenshannon(aucs/aucs.sum(), labels/labels.sum())

    for cidx, regulon_name in enumerate(regulons):
        for ridx, cell_type in enumerate(cell_types):
            rss_values[ridx, cidx] = rss(auc_mtx[regulon_name], (cell_type_series == cell_type).astype(int))

    return pd.DataFrame(data=rss_values, index=cell_types, columns=regulons)


rss=regulon_specificity_scores(auc_mtx, anno['cluster']).T
rss.to_csv('data/scenic/HT280_regulon_specificity_scores.csv')
