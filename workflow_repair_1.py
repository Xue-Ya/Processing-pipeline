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
# print(sys.path)
from Repair import *
rcParams['pdf.fonttype'] = 42

# import scvi
# print(scvi.__version__)

# RNA_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Assembled10Domains.h5ad'
celltype_key = 'new_celltype'

sc_meta_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/sc_meta.csv'
st_meta_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/st_meta.csv'

sc_gene_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/sc_gene.csv'
st_gene_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/st_gene.csv'

output_path = 'FigureData/Dataset1/stage2/Result_Repair/'
conv_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/stage2/result_deconv/Tangram_result.csv'

if not os.path.exists(output_path):
    os.makedirs(output_path)

test = Repair(RNA_file = sc_gene_path, Spatial_file = st_gene_path, Spitial_meta = st_meta_path, celltype_file = sc_meta_path, output_path = output_path,conv_path = conv_path,celltype_key = celltype_key)
Methods = ['spatalk', 'sprout']
Result = test.Repair(Methods)