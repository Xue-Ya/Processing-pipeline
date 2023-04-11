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

import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

root_dir = os.path.join(os.getcwd(),'')
if root_dir not in sys.path:
    sys.path.append(root_dir)
logger.info(sys.path)
from DeconvolutionSpot import *
from Repair import *
rcParams['pdf.fonttype'] = 42

# dataset = "bulk_dataset1"
dataset = "Dataset3"
# dataset = "Dataset1"
overwrite = False
# RNA_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/res_filtered.h5ad'
# RNA_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/res_healthy_sub.h5ad'
# RNA_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/res_cancer_sub.h5ad'
RNA_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/res_healthy.h5ad'
# RNA_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/res_cancer.h5ad'
# RNA_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/res_subsample.h5ad'

st = sys.argv[1]

    # celltype_key = sys.argv[3]
    # output_path = sys.argv[4]
    # Methods = sys.argv[5]
    
# if dataset == "'+dataset + '":   
work_dict = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/'+ dataset + '/apart/'
        
output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/'+dataset + '/' + st +'/Result_Deconv/'

data_save = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC'
# else:
#     work_dict = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/'+dataset + '/filtered/'
            
#     output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/'+dataset + '/'+ dorner +'/Result_Deconv/'
    
if not os.path.exists(data_save):
    os.makedirs(data_save)

if not os.path.exists(data_save+ '/Spatalk/'):
    os.makedirs(data_save+ '/Spatalk/')

if not os.path.exists(output_path):
    os.makedirs(output_path)
    
if not os.path.exists(work_dict + st + '/Spatalk/'):
    os.makedirs(work_dict + st + '/Spatalk/')


# RNA_h5ad = '/home/xuezhengyang/data6/02-deconv_1/Script/Assembled10Domains.h5ad'
Spatial_h5ad = work_dict + 'st/' + st + '.h5ad'

# st_adata =  '/home/xuezhengyang/data6/02-deconv_1/Script/Data/'+dataset + '/apart/st/st_12.h5ad'
count_f =  work_dict + st + '/Spatalk/spatalk_st_gene.csv'
meta_f = work_dict + st + '/Spatalk/spatalk_st_corr.csv'
# sc_adata = sc.read('/home/xuezhengyang/data6/02-deconv_1/Script/Data/'+dataset + '/apart/sc/stage2.h5ad')
# sc_count = data_save + '/Spatalk/spatalk_sc_gene_sub.csv'
# sc_celltypes = data_save + '/Spatalk/spatalk_sc_celltype_sub.csv'

sc_count = data_save + '/Spatalk/spatalk_sc_gene_healthy.csv'
sc_celltypes = data_save + '/Spatalk/spatalk_sc_celltype_healthy.csv'

# sc_count = data_save + '/Spatalk/spatalk_sc_gene_cancer.csv'
# sc_celltypes = data_save + '/Spatalk/spatalk_sc_celltype_cancer.csv'

# sc_gene_path =  data_save + '/sc_gene_sub.csv'
# sc_meta_path = data_save + '/sc_meta_sub.csv'

sc_gene_path =  data_save + '/sc_gene_healthy.csv'
sc_meta_path = data_save + '/sc_meta_healthy.csv'

# sc_gene_path =  data_save + '/sc_gene_cancer.csv'
# sc_meta_path = data_save + '/sc_meta_cancer.csv'

st_corr_path = work_dict + st + '/st_corr.csv'
st_gene_path = work_dict + st + '/st_gene.csv'

logger.info(st + ' spatalk data prepare begins')

os.system('python /home/xuezhengyang/data6/02-deconv_1/Script/data_prepare/Script1.py ' + Spatial_h5ad + ' '+ count_f + ' '+ meta_f + ' '+ RNA_h5ad + ' '+ sc_count + ' '+ sc_celltypes )

logger.info(st + ' spatalk data prepare is done')

logger.info( st + ' other data prepare begins')
os.system('python /home/xuezhengyang/data6/02-deconv_1/Script/data_prepare/Script2.py ' + RNA_h5ad + ' '+ Spatial_h5ad + ' '+ sc_meta_path + ' '+ st_corr_path + ' '+ sc_gene_path + ' '+ st_gene_path )
logger.info(st + ' other data prepare is done')

# celltype_key = 'celltype_final'
celltype_key = 'cell_type'
# celltype_key = 'new_celltype'
Methods = []
test = Deconvolutions(RNA_h5ad = RNA_h5ad, Spatial_h5ad = Spatial_h5ad, celltype_key = celltype_key, output_path = output_path,work_dir= work_dict + st,sc_dir= data_save)

	
# if ((not os.path.exists(output_path + 'DestVI_result.csv')) or (overwrite)):
#     Methods.append('DestVI')
if  ((not os.path.exists(output_path + 'Tangram_result.csv')) or (overwrite)):
    Methods.append('Tangram')
# if  ((not os.path.exists(output_path + 'SpaTalk_result.csv')) or (overwrite)):
#     Methods.append('SpaTalk')
# if  ((not os.path.exists(output_path + 'RCTD_result.csv')) or (overwrite)):
#     Methods.append('RCTD')
# if  ((not os.path.exists(output_path + 'Seurat_result.csv')) or (overwrite)):
#     Methods.append('Seurat')

# if  ((not os.path.exists(output_path + 'SPOTlight_result.csv')) or (overwrite)):
#     Methods.append('SPOTlight')
# if  ((not os.path.exists(output_path + 'deconvSeq_result.csv')) or (overwrite)):
#     Methods.append('deconvSeq')
# if  ((not os.path.exists(output_path + 'stereoscope_result.csv')) or (overwrite)):
#     Methods.append('stereoscope')
# if  ((not os.path.exists(output_path + 'cell2location_result.csv')) or (overwrite)):
#     Methods.append('cell2location')
# Tangram',
# Result = test.Dencon(Methods)

logger.info('Deconv is Done')



# 'SpaTalk','RCTD','Seurat','SPOTlight','deconvSeq','stereoscope','cell2location'
spa_methods = ['RCTD','Seurat','SPOTlight','deconvSeq','stereoscope','cell2location']
DMethods = ['Tangram', 'Seurat', 'Seurat']
# 'SpaTalk',
# 'stereoscope',
# 'Tangram','SpaTalk','RCTD','Seurat','SPOTlight',

for method in DMethods:    
    conv_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/'+dataset + '/' + st +'/Result_Deconv/' + method + '_result.csv'
    not_exist =  (not os.path.exists(conv_path))
    
    if (method == "SpaTalk") and not_exist:
        conv = 1
        Methods = ['spatalk']
    elif method == "RCTD" and not_exist:
        conv = 2
        Methods = ['spatalk']
    elif method == "Seurat" and not_exist:
        Methods = ['spatalk']
        conv = 3
    elif method == "SPOTlight" and not_exist:
        Methods = ['spatalk']
        conv = 4
    elif method == "deconvSeq" and not_exist:
        Methods = ['spatalk']
        conv = 5
    elif method == "stereoscope" and not_exist:
        Methods = ['spatalk']
        conv = 6
    elif method == "cell2location" and not_exist:
        Methods = ['spatalk']
        conv = 7
    else:
        conv = 0
        if method in spa_methods:
            Methods = ['sprout']
        else:
            Methods = ['sprout','spatalk']
    
    
    #     if (method == "SpaTalk") and not_exist:
    #     conv = 1
    #     Methods = ['spatalk']
    # elif method == "RCTD" and not_exist:
    #     conv = 2
    #     Methods = ['spatalk']
    # elif method == "Seurat" and not_exist:
    #     Methods = ['spatalk']
    #     conv = 3
    # elif method == "SPOTlight" and not_exist:
    #     Methods = ['spatalk']
    #     conv = 4
    # elif method == "deconvSeq" and not_exist:
    #     Methods = ['spatalk']
    #     conv = 5
    # elif method == "stereoscope" and not_exist:
    #     Methods = ['spatalk']
    #     conv = 6
    # elif method == "cell2location" and not_exist:
    #     Methods = ['spatalk']
    #     conv = 7

    # Methods = ['sprout','spatalk']
    
    output_path = 'FigureData/'+ dataset + '/' + st +'/Result_Repair/' + method

    # output_path = 'FigureData/'+dataset + '/' + st +'/Result_Repair/DestVI/Tangram'
    # conv_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/'+dataset + '/'+ rna +'_' + st +'/Result_Deconv/Tangram_result.csv'

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if not os.path.exists(output_path + '/sprout'):
        os.makedirs(output_path + '/sprout')

    if not os.path.exists(output_path + '/spatalk'):
        os.makedirs(output_path + '/spatalk')

    output_path = output_path + '/'
    test = Repair(sc_dir= data_save,work_dir= work_dict + st , output_path = output_path,conv_path = conv_path,celltype_key = celltype_key, conv=conv)
    
    # 'sprout',
    Result = test.Repair(Methods)
    
    print(method + ' recon done')