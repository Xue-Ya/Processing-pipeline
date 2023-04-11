import scanpy as sc
import pandas as pd
import anndata as ad

scST_path = ''
st_path = ''
sc_path = ''
corr_path = ''
sc_cell_type = 'cell_type'
st_Cell_type = 'annotation'

scst = pd.read_csv(scST_path,index_col=0)

st = sc.read(st_path)

sc_adata = sc.read(sc_path)

scst.index = scst['index']
scst = scst.drop(['index'],axis=1)

celltype = set(sc_adata.obs[sc_cell_type].index)

cell = scst.index
index = pd.DataFrame(index=scst.index)

index[sc_cell_type] = sc_adata.obs[sc_cell_type].loc[index.index]

scst_adata = ad.AnnData(scst,obs = index, dtype='int32')

corr = pd.read_csv(corr_path,sep='\t',index_col=0)

corr = corr.values[:]

tmp = corr[:,1].copy()
corr[:,1] = corr[:,0]
corr[:,0] = tmp

scst_adata.obsm['spatial'] = corr


####scST gene map 
sc.pl.spatial(scst_adata, color= 'new_celltype' ,scale_factor = 0.1,spot_size = 0.5)


####ST gene map
sc.pl.spatial(st,color='annotation')

