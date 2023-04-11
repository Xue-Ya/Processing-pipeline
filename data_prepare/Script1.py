import scanpy as sc
import sys
import pandas as pd
from os.path import exists

# st_path = sys.argv[]
# sc_path = sys.argv[]
# count_f = sys.argv[]
# meta_f = sys.argv[]


# output_path = sys.argv[]
# output_path = ''
overwrite = False

if ((not exists(sys.argv[2])) or (not exists(sys.argv[3])) or overwrite):
    st_adata =  sc.read(sys.argv[1])
    # st_adata = sc.read('/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/st/st_12.h5ad')
    cluster_key = 'annotation'
    pos_key = ['array_row','array_col']
    is_human = True
    count_f = sys.argv[2] 
    # count_f = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/Spatalk/spatalk_Dataset1_st_gene.csv'
    meta_f = sys.argv[3] 
    # meta_f = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/Spatalk/spatalk_st_meta.csv'
    df = st_adata.to_df().T
    # df.index = "C"+df.index
    
    df.to_csv(count_f)
    meta = st_adata.obs[pos_key+[cluster_key]].reset_index()
    meta.columns = ['spot', 'x', 'y', 'celltype']
    # meta.spot = "C"+meta.spot
    if not pd.api.types.is_string_dtype(meta.celltype.dtype):
        # meta.celltype = "T"+meta.celltype.astype('str')
        meta.celltype = meta.celltype.astype('str')
    meta.to_csv(meta_f)
    species = 'Human' if is_human else 'Mouse'

if (not exists(sys.argv[5])):
    sc_adata = sc.read(sys.argv[4])
    # sc_adata = sc.read('/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/sc/stage2.h5ad')
    sc_df = sc_adata.to_df().T
    sc_count = sys.argv[5] 
    # sc_count = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/Spatalk/spatalk_Dataset1_sc_gene.csv'
    # sc_df.index = "C"+sc_df.index
    sc_df.to_csv(sc_count)
    
if (not exists(sys.argv[6])):
    sc_celltypes = sys.argv[6] 
    # sc_celltypes = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/Spatalk/sc_celltype.csv'
    # sc_meta = sc_adata.obs['new_celltype'].reset_index()
    sc_meta = sc_adata.obs['cell_type'].reset_index()
    sc_meta.columns = ['cell', 'celltype']
    # sc_meta.gene = "C"+sc_meta.gene
    if not pd.api.types.is_string_dtype(sc_meta.celltype.dtype):
        # sc_meta.celltype = "T"+sc_meta.celltype.astype('str')
        sc_meta.celltype = sc_meta.celltype.astype('str')
    sc_meta.to_csv(sc_celltypes)