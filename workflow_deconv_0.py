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
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr,ttest_ind,mannwhitneyu
from sklearn.metrics import mean_squared_error

root_dir = os.path.join(os.getcwd(),'')
if root_dir not in sys.path:
    sys.path.append(root_dir)
print(sys.path)
from DeconvolutionSpot import *
rcParams['pdf.fonttype'] = 42

# RNA_h5ad = sys.argv[1]
# Spatial_h5ad = sys.argv[2]
# celltype_key = sys.argv[3]
# output_path = sys.argv[4]
# Methods = sys.argv[5]

RNA_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/sc/stage2.h5ad'
# RNA_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Assembled10Domains.h5ad'
Spatial_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/st/st_12.h5ad'
# celltype_key = 'celltype_final'
celltype_key = 'new_celltype'
output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/stage2/result_deconv/'
if not os.path.exists(output_path):
    os.makedirs(output_path)

test = Deconvolutions(RNA_h5ad = RNA_h5ad, Spatial_h5ad = Spatial_h5ad, celltype_key = celltype_key, output_path = output_path)
Methods = ['Tangram','DestVI']
Result = test.Dencon(Methods)