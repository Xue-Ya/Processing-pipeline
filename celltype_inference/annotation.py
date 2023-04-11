#!/usr/bin/env /home/xuezhengyang/.conda/envs/FYP/bin/python
import pandas as pd
import scanpy as sc
import anndata as ad

# rna = "stage2"
# st = "st_12"

# rna="stage2"
# rna="stage3"
# regions=["st_12","st_12_2","st_14","st_14_2"]
# regions = ['st_12']
# ,"st_12_2","st_14","st_14_2"
samples = ["sample1","sample2","sample3","sample4","sample5","sample6",
           "sample7","sample8","sample9","sample10","sample11"]
# regions=[ "_TD2_","_TD5_","_TD6_","_TD8_","_TD3_"]
# regions = ['st_20_3']
deconv = 'Tangram'
# res = 'spatalk'
dataset = "Dataset3"
# dataset = "dataset1"

for i in samples:
    st_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/'+ dataset +'/apart/st/'+ i +'.h5ad'
    
    st = sc.read(st_path)
    
    conv_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/'+dataset + '/' + i +'/Result_Deconv/' + deconv + '_result.csv'
    df = pd.read_csv(conv_path,index_col=0)
    
    annotation = df.idxmax(axis=1)
    
    st.obs['annotation'] = annotation
    
    st.write_h5ad(st_path)