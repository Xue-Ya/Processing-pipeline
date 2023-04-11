import pandas as pd
import scanpy as sc


def draw_st_spot(st,save_path,dataset,celltype_key,pois_key):
    
    data_dir = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/'+dataset+'/apart/'
    # save_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Figure/' + rna+'_'+st + '/st_deconv'
    
    st_adata = sc.read(data_dir + 'st/' + st + '.h5ad')

    # row = pd.DataFrame(st.obs[pois_key[0]])
    # col = pd.DataFrame(st.obs[pois_key[1]])

    # corr = pd.concat([row,col], axis=1)
    corr = pd.read_csv(data_dir + st +'/st_corr.csv', index_col=0)
    
    meta = pd.DataFrame(st_adata.obs[celltype_key])

    res = pd.concat([corr,meta], axis=1)
    
    res.index = res.index.str.replace('-','_')
    res.index.name='id'
    
    res.to_csv(save_path + 'st_spot.csv')