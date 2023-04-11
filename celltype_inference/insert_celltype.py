#!/usr/bin/env /home/xuezhengyang/.conda/envs/FYP/bin/python

import sys
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import pandas as pd
import scanpy as sc
import numpy as np
import anndata as ad

path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/'

#SC
# samples = ['GSM5699777','GSM5699778','GSM5699779','GSM5699781','GSM5699782','GSM5699784']

#ST
samples = ['GSM5702473', 'GSM5702474','GSM5702475','GSM5702476','GSM5702477','GSM5702478']

dorners = ['_TD1_','_TD2_','_TD3_','_TD5_','_TD6_','_TD8_']

for i in range(0,len(samples)):
    dorner = dorners[i]
    
    # adata = sc.read(path +'sc/' + dorner + '.h5ad')
    adata = sc.read(path +'st/' + dorner + '.h5ad')
    
    # cell_type = pd.read_csv(path +'filtered/sc/sc' + dorner + '.csv',index_col=0)
    cell_type = pd.read_csv(path +'filtered/st/st' + dorner + '.csv',index_col=0)
    
    # adata.obs['new_celltype'] = cell_type
    adata.obs['annotation'] = cell_type
    
    # adata.write(path +'filtered/sc/sc' + dorner + '.h5ad')
    adata.write(path +'filtered/st/st' + dorner + '.h5ad')
    
    