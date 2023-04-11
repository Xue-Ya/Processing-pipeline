import pandas as pd
import scanpy as sc
import anndata as ad

def scST_exp(st,repair_path,dataset,save_path):
    
    sc_cell_type = 'new_celltype'
    st_Cell_type = 'annotation'
    data_dir = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/'+dataset+'/apart/'
    
    
    scst = pd.read_csv(repair_path+ 'sprout/sc_agg_exp.tsv',index_col=0)
    
    scst.index = scst['index']
    scst = scst.drop(['index'],axis=1)
    
    cell = scst.index
    index = pd.DataFrame(index=scst.index)

    scst_adata = ad.AnnData(scst,obs = index, dtype='int32')
    
    exp = scst_adata.to_df()
    
    exp.index.name='id'
    
    exp.to_csv(save_path + 'scST_exp.csv')
    