import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import os
from scipy.stats import pearsonr
import script4

sc_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/sc/stage2.h5ad'
st_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/st/st_12.h5ad'
bulk_gct = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset0/gene_tpm_2017-06-05_v8_lung.gct'
output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2_st_12/bulk'

if not os.path.exists(output_path):
    os.makedirs(output_path)

def generate_data(sc_path, st_path, bulk_gct,output_path, sc_Cell_type = 'new_celltype', st_Cell_type = 'annotation', pois_key = ['array_row','array_col']):
    sc_data = sc.read(sc_path)
    st_data = sc.read(st_path)
    
    sc_data = script4.filter(sc_data, sc_Cell_type)
    st_data = script4.filter(st_data, st_Cell_type)
    # invalid_groups_selected = sc_data.obs[sc_Cell_type].value_counts().loc[lambda x: x < 2].index
    # cell = pd.DataFrame(sc_data.obs[sc_Cell_type])

    # cell_filter = cell.loc[cell[sc_Cell_type].isin(invalid_groups_selected)]

    # sc_data = sc_data[~sc_data.obs_names.isin(cell_filter.index)]
    
    # new_cell = cell.loc[~cell[sc_Cell_type].isin(invalid_groups_selected)]
    # sc_data.obs[sc_Cell_type] = pd.Series(new_cell[sc_Cell_type].values, index=new_cell.index)
    
    # cell_filter = pd.Series(new_cell[st_Cell_type], index=new_cell.index,dtype=str)
    # cell_filter = pd.Series(cell_filter, index=cell_filter.index,dtype='category')

    # sc_data.obs[st_Cell_type] = cell_filter
    
    #Generate sc_data.csv
    sc_data_X = sc_data.to_df().T
    
    #Generate st_data.csv
    
    st_data_X =  st_data.to_df().T   
    
    #Generate sc_meta.csv
    sc_meta = pd.DataFrame(sc_data.obs[sc_Cell_type])
    # sc_meta.to_csv(output_path + '/sc_meta.csv')
    # sc_meta = pd.read_csv(output_path + '/sc_meta.csv')
    list = []
    list.append("Cell_type")
    sc_meta.columns = list
    
    print(sc_meta)
    sc_meta['Cell'] = sc_meta.index

    order = ['Cell', "Cell_type"]
    sc_meta = sc_meta[order]
    sc_meta.to_csv(output_path + '/sc_meta.csv')
    
    
    #Generate st_meta.csv
    st_meta = st_data.obs[[st_Cell_type]+pois_key]
    list = ['Cell_type','xcoord','ycoord']
    st_meta.columns = list
    st_meta['Cell'] = st_meta.index

    order = ['Cell', "Cell_type",'xcoord','ycoord']
    st_meta = st_meta[order]
    st_meta.to_csv(output_path + '/st_meta.csv')
    
    #Generate bulk.csv
    bulk_data = pd.read_csv(bulk_gct, sep='\t',skiprows=2,index_col=2)
    bulk_data = bulk_data.drop("Name",axis=1)
    bulk_data = bulk_data.drop("id",axis=1)
    bulk_data = bulk_data.groupby(level=0).last()
    
    
    sc_data_X.to_csv(output_path + '/sc_data.csv')
    # sc_meta = pd.read_csv(output_path + '/sc_meta.csv')
    st_data_X.to_csv(output_path + '/st_data.csv')
    # st_meta = pd.read_csv(output_path + '/st_meta.csv')
    bulk_data.to_csv(output_path + '/bulk.csv')
    bulk_data.iloc[:,0].to_csv(output_path + '/bulk_demo.csv')
    
generate_data(sc_path, st_path,bulk_gct, output_path)