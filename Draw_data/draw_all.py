#!/usr/bin/env /home/xuezhengyang/.conda/envs/FYP/bin/python
import sys
sys.path.append("/home/xuezhengyang/data6/02-deconv_1/Script/Draw_data")

import decon
import st_spot
import scST_spot
import marker_genes
import scST_exp
import pandas as pd
import os,sys
import scanpy as sc
import shutil
import spatalk

# regions=["st_12","st_12_2","st_14","st_14_2"]
# ,"st_12_2","st_14","st_14_2"
regions=["st_12","st_12_2","st_14","st_14_2","st_16","st_16_2","st_17","st_17_2"]
dataset = 'Dataset1'
deconvs = ['Tangram','Seurat']
# ,"st_20_3"

st_celltype_key = 'annotation'
pois_key = ['array_row','array_col']
dataset_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/'+dataset+'/apart/'

for st in regions:
    print(st)
    output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Figure/' + st

    st_spots = output_path + '/st_spot'
    st_deconv = output_path + '/st_deconv'
    
    for deconv in  deconvs:
        scst = output_path + '/scst/' + deconv 
        marker_gene = output_path + '/marker_gene/' + deconv 
        scst_gene = output_path + '/scst_gene/' + deconv 
        marker_spatalk =  output_path + '/marker_gene/'+ deconv + '/spatalk'
        
        if not os.path.exists(scst):
            os.makedirs(scst)
            
        if not os.path.exists(marker_gene):
            os.makedirs(marker_gene)
            
        if not os.path.exists(scst_gene):
            os.makedirs(scst_gene)
        
        if not os.path.exists(marker_spatalk):
            os.makedirs(marker_spatalk)
    
    if not os.path.exists(st_spots):
        os.makedirs(st_spots)
    
    if not os.path.exists(st_deconv):
        os.makedirs(st_deconv)
        


#Generate Deconv file
for st in regions:
    data_dir = dataset_path + st + '/'
    res_dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/'+dataset+'/' + st
    decon_path = res_dir + '/Result_Deconv/'
    save_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Figure/' + st + '/st_deconv/'
    methods = ['Tangram_result.csv','Seurat_result.csv']

    # decon.draw_deconv(data_dir,save_path,decon_path,methods)

#Generate st_spot file
for st in regions:
    data_dir = dataset_path + st + '/'
    save_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Figure/' + st + '/st_spot/'
    
    st_spot.draw_st_spot(st,save_path,dataset,st_celltype_key,pois_key)

#Generate scST_spot sprout file
for st in regions:
    res_dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/'+dataset+'/' + st
    repair_path = res_dir + '/Result_Repair/'
    save_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Figure/' + st + '/scst/'
    
    # scST_spot.draw_scst_spot(st,save_path,repair_path)

#Generate marker gene file
for st in regions:
    res_dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/'+dataset+'/' + st
    repair_path = res_dir + '/Result_Repair/'
    save_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Figure/' + st + '/marker_gene/'
    
    # marker_genes.draw_marker_genes(st,repair_path=repair_path,dataset = dataset, save_path= save_path)
    
for st in regions:
    res_dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/'+dataset+'/' + st
    repair_path = res_dir + '/Result_Repair/'
    save_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Figure/' + st + '/scst_gene/'
    
    # scST_exp.scST_exp(st,repair_path=repair_path,dataset = dataset,save_path= save_path)
    
for st in regions:
    save_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Figure/' + st
    
    
    for deconv in deconvs:
        spatalk.spatalk(st,deconv,dataset,save_path)
