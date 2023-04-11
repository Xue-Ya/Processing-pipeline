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
# import DeconvolutionSpot
rcParams['pdf.fonttype'] = 42
plt.style.use('default')

# path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/'
# sample = 'GSM5699784'
# dorner = '_TD8_'

# barcode = pd.read_csv(path +'sc/'+ sample + dorner + 'barcodes.tsv',sep='\t')
# features = pd.read_csv(path +'sc/'+ sample + dorner + 'features.tsv',sep='\t')
# matrix = sc.read(path +'sc/'+ sample + dorner + 'matrix.mtx')

def col(barcode):
    barcode = barcode.T

    barcode = barcode.reset_index()
    barcode = barcode.T
    barcode = barcode.reset_index(drop=True)
    return barcode


path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/'

# samples = ['GSM5699777','GSM5699778','GSM5699779','GSM5699781','GSM5699782','GSM5699784']
samples = ['GSM5702473', 'GSM5702474','GSM5702475','GSM5702476','GSM5702477','GSM5702478']
dorners = ['_TD1_','_TD2_','_TD3_','_TD5_','_TD6_','_TD8_']

# for i in range(0,len(samples)):
#     sample = samples[i]
#     dorner = dorners[i]

#     # barcode = pd.read_csv(path +'sc/'+ sample + dorner + 'barcodes.tsv',sep='\t')
#     # features = pd.read_csv(path +'sc/'+ sample + dorner + 'features.tsv',sep='\t')
#     # matrix = sc.read(path +'sc/'+ sample + dorner + 'matrix.mtx')
    
#     barcode = pd.read_csv(path +'st/'+ sample + dorner + 'barcodes.tsv',sep='\t')
#     features = pd.read_csv(path +'st/'+ sample + dorner + 'features.tsv',sep='\t')
#     matrix = sc.read(path +'st/'+ sample + dorner + 'matrix.mtx')

#     tissue_pos = pd.read_csv(path +'st/'+ sample + dorner + 'tissue_positions_list.csv')
    
#     pos_filter = tissue_pos.loc[tissue_pos.iloc[:,0] == 1]

#     x = pos_filter.iloc[:,3]
#     y = pos_filter.iloc[:,4]

#     # matrix = matrix[0:-1,0:-1]

#     barcode = col(barcode)
#     features = col(features)

#     matrix.var_names = barcode.values[:,0]
#     matrix.obs_names = features.values[:,1] 

#     matrix = matrix.T

#     matrix.obs['array_row'] = x
#     matrix.obs['array_col'] = y

#     matrix.obs['batch'] = dorner
    
#     matrix.write(path +'st/' + dorner + '.h5ad')


adata_new = sc.read(path +'st/_TD1_.h5ad')
adata_new.var_names_make_unique()
for dorner in dorners[1:]:
    sc_new = sc.read(path +'st/' + dorner + '.h5ad')
    sc_new.var_names_make_unique()
    adata_new = ad.concat([adata_new,sc_new],uns_merge = "first")

adata_new.obs['batch'].unique()
sc.pp.combat(adata_new,inplace = True)

adata_new.write(path +'st/all.h5ad')