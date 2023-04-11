import scanpy as sc
import sys
import pandas as pd
from os.path import exists


# sc_adata = sc.read('/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/sc/stage2.h5ad')
# st_adata = sc.read('/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/st/st_12.h5ad')

# sc_meta_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/sc_meta.csv'
# st_meta_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/st_meta.csv'

# sc_gene_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/sc_gene.csv'
# st_gene_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/stage2/st_gene.csv'

# output_path = 'FigureData/Dataset1/Result_Repair/stage2/'


st_adata = sc.read(sys.argv[2])
sc_meta_path = sys.argv[3]
st_corr_path = sys.argv[4]
sc_gene_path = sys.argv[5]
st_gene_path = sys.argv[6]

if ((not exists(sc_gene_path))|(not exists(sc_meta_path))):
    sc_adata = sc.read(sys.argv[1])
    sc_df = sc_adata.to_df()
    sc_df.to_csv(sc_gene_path)
    # sc_meta = sc_adata.obs['new_celltype']
    sc_meta = sc_adata.obs['cell_type']
    sc_meta.to_csv(sc_meta_path)

if ((not exists(st_gene_path))|(not exists(st_corr_path))):
    df = st_adata.to_df()
    df.to_csv(st_gene_path, index_label= False)
    row = pd.DataFrame(st_adata.obs['array_row'])
    col = pd.DataFrame(st_adata.obs['array_col'])
    cell_type = pd.DataFrame(st_adata.obs['annotation'])
    corr = pd.concat([col,row], axis=1)
    corr.to_csv(st_corr_path)






