import sys
import os

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc


PLOT_DIR = './figures/'


def load_BD_csvs(sample_list):
    adata_list = []
    for sample in sample_list:
        filedir = f'./data/No3_1_{sample}_RSEC_MolsPerCell.csv'
        adata = load_csv(filedir, sample=sample, comment='#')
        adata_list.append(adata)
    adata = ad.concat(adata_list)
    return adata


def load_csv(filedir, sample, comment=''):
    df = pd.read_table(filedir, sep=',', comment=comment, index_col=0)
    adata = sc.AnnData(X=df)
    adata.obs['dataset'] = sample
    return adata
