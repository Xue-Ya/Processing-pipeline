#!/usr/bin/env /home/xuezhengyang/.conda/envs/FYP/bin/python
import scanpy as sc
import pandas as pd

st  = True

def col(barcode):
    barcode = barcode.T

    barcode = barcode.reset_index()
    barcode = barcode.T
    barcode = barcode.reset_index(drop=True)
    return barcode

path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/'


if st:
    samples = ['GSM5702473','GSM5702474','GSM5702475','GSM5702476','GSM5702477','GSM5702478']
    dorners = ['_TD1_','_TD2_','_TD3_','_TD5_','_TD6_','_TD8_']
    for i in range(0,len(samples)):
        sample = samples[i]
        dorner = dorners[i]
        
        barcode = pd.read_csv(path +'st/'+ sample + dorner + 'barcodes.tsv',sep='\t')
        features = pd.read_csv(path +'st/'+ sample + dorner + 'features.tsv',sep='\t')
        matrix = sc.read(path +'st/'+ sample + dorner + 'matrix.mtx')
        tissue_pos = pd.read_csv(path +'st/'+ sample + dorner + 'tissue_positions_list.csv',index_col=0)

        # matrix = matrix[0:-1,0:-1]

        barcode = col(barcode)
        features = col(features)

        matrix.var_names = barcode.values[:,0]
        matrix.obs_names = features.values[:,1]

        pos_filter = tissue_pos.loc[tissue_pos.iloc[:,0] == 1]

        x = pos_filter.iloc[:,3]
        y = pos_filter.iloc[:,4]

        matrix = matrix.T
        
        matrix = matrix[list(pos_filter.index)]
        
        matrix.obs['array_row'] = x
        matrix.obs['array_col'] = y
        
        matrix.obs['batch'] = dorner

        matrix.write(path +'st/' + dorner + '.h5ad')
else:
    samples = ['GSM5699777','GSM5699778','GSM5699779','GSM5699781','GSM5699782','GSM5699784']
    dorners = ['_TD1_','_TD2_','_TD3_','_TD5_','_TD6_','_TD8_']
    
    for i in range(0,len(samples)):
        sample = samples[i]
        dorner = dorners[i]
        
        barcode = pd.read_csv(path +'sc/'+ sample + dorner + 'barcodes.tsv',sep='\t')
        features = pd.read_csv(path +'sc/'+ sample + dorner + 'features.tsv',sep='\t')
        matrix = sc.read(path +'sc/'+ sample + dorner + 'matrix.mtx')

        # matrix = matrix[0:-1,0:-1]

        barcode = col(barcode)
        features = col(features)

        matrix.var_names = barcode.values[:,0]
        matrix.obs_names = features.values[:,1]

        matrix.write(path +'sc/'+ sample + dorner + '.h5ad')