import pandas as pd
import scanpy as sc

# rna="stage2"
# # rna="stage3"
# regions=("st_12","st_12_2","st_14","st_14_2")
# regions=["st_16","st_16_2","st_17","st_17_2","st_19","st_20_1","st_20_2","st_20_3"]

def draw_deconv(data_dir,save_path,decon_path,methods):


    corr = pd.read_csv(data_dir + 'st_corr.csv', index_col=0)
    
    corr.index = corr.index.str.replace('-','_')
    
    for method in methods:
        data = pd.read_csv(decon_path + method,index_col=0)
        data.index = data.index.str.replace('-','_')
        res = pd.concat([corr,data],axis= 1)
        res.rename(columns={'array_row':'Spot_xcoord'},inplace=True)
        res.rename(columns={'array_col':'Spot_ycoord'},inplace=True)
        res.index.name='id'
        res.to_csv(save_path + method)
