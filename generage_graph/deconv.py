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
# from scipy.spatial.distance import jensenshannon
# from scipy.stats import pearsonr,ttest_ind,mannwhitneyu
# from sklearn.metrics import mean_squared_error
# import DeconvolutionSpot
# import os,sys
root_dir = os.path.join(os.getcwd(),'')
if root_dir not in sys.path:
    sys.path.append(root_dir)
    
rcParams['pdf.fonttype'] = 42
plt.style.use('default')
import anndata as ad
import squidpy as sq

dataset = 'Dataset1'
dataset = 'bulk_dataset1'

# sts = ["st_12_2","st_14","st_14_2","st_12","st_16","st_16_2","st_17","st_17_2","st_19","st_20_1","st_20_2","st_20_3"]
sts = ["_TD1_", "_TD2_", "_TD5_", "_TD6_", "_TD8_"]

out = '/home/xuezhengyang/data6/02-deconv_1/Script/figures/show/st_decon/' + dataset + '/scst/'

if not os.path.exists(out):
    os.makedirs(out)
    
for st in sts:
    deconv = pd.read_csv('/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/' + dataset + '/' + st + '/Result_Deconv/Tangram_result.csv',index_col=0)
    adata_st = sc.read('/home/xuezhengyang/data6/02-deconv_1/Script/Data/' + dataset + '/apart/st/' + st + '.h5ad')
    
    adata_st.obs = pd.concat([adata_st.obs, deconv], axis=1)
    
    sc.pl.spatial(
    adata_st,
    # color = 'annotation'
    color=deconv.columns,
    save = '/st_decon/' + dataset + '/scst/' +st + '.png',
    show=False
    # , "CD4-positive, alpha-beta T cell", "pericyte", "ionocyte", "type II pneumocyte", "B cell", "type I pneumocyte"
)
    