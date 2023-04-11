#!/usr/bin/env /home/xuezhengyang/.conda/envs/FYP/bin/python
import scanpy as sc
import pandas as pd
import anndata as ad

out_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1'
work_dir = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1'
save = '/home/xuezhengyang/data6/02-deconv_1/Script/Figure'
save = save + '/stage2_st_12/marker_gene'
# /home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/stage2_st_12/Result_Repair/sprout/sc_agg_exp.tsv
scST_path = out_path + '/stage2_st_12/Result_Repair/sprout/sc_agg_exp.tsv'
st_path = work_dir + '/apart/st/st_12.h5ad'
sc_path = work_dir + '/apart/sc/stage2.h5ad'
corr_path = out_path + '/stage2_st_12/Result_Repair/spot_raw_coord_best.tsv'
# /home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/stage2_st_12/Result_Repair/spot_raw_coord_best.tsv
sc_cell_type = 'new_celltype'
st_Cell_type = 'annotation'

scst = pd.read_csv(scST_path, index_col=0)

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

from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


############ scST marker gene export

cell = pd.DataFrame(scst_adata.obs[sc_cell_type])

cell_filter = pd.Series(cell[sc_cell_type], index=cell.index,dtype=str)
cell_filter = pd.Series(cell_filter, index=cell_filter.index,dtype='category')

scst_adata_new = scst_adata.copy()

scst_adata_new.obs[sc_cell_type] = cell_filter

sc.pp.log1p(scst_adata_new)

sc.tl.rank_genes_groups(scst_adata_new, sc_cell_type,  method='wilcoxon')

sc.pl.rank_genes_groups(scst_adata_new,sharey=False,n_genes=10,show = False)

scst_marker_gene = pd.DataFrame(scst_adata_new.uns['rank_genes_groups']['names']) .head(5)

result = scst_adata_new.uns['rank_genes_groups']
groups = result['names'].dtype.names
scst_p = pd.DataFrame(
    {group + ' '+ key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)


############ ST marker gene export

st_Cell_type = 'annotation'

st_new = st.copy()

invalid_groups_selected = st_new.obs[st_Cell_type].value_counts().loc[lambda x: x < 2].index
cell = pd.DataFrame(st_new.obs[st_Cell_type])

cell_filter = cell.loc[cell[st_Cell_type].isin(invalid_groups_selected)]

st_new = st_new[~st_new.obs_names.isin(cell_filter.index)]

new_cell = cell.loc[~cell[st_Cell_type].isin(invalid_groups_selected)]

st_new.obs[st_Cell_type] = pd.Series(new_cell[st_Cell_type].values, index=new_cell.index)

cell_filter = pd.Series(new_cell[st_Cell_type], index=new_cell.index,dtype=str)
cell_filter = pd.Series(cell_filter, index=cell_filter.index,dtype='category')

st_new.obs[st_Cell_type] = cell_filter

sc.pp.log1p(st_new)

sc.tl.rank_genes_groups(st_new, st_Cell_type,  method='wilcoxon')

sc.pl.rank_genes_groups(st_new,sharey=False,n_genes=10,show = False)

st_marker_gene = pd.DataFrame(st_new.uns['rank_genes_groups']['names']).head(5)

result = st_new.uns['rank_genes_groups']
groups = result['names'].dtype.names
st_p = pd.DataFrame(
    {group + ' '+ key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)

########## To csv

import os

if not os.path.exists(save):
    os.makedirs(save)
    
scst_marker_gene.to_csv(save + 'scst_marker.csv')

st_marker_gene.to_csv(save + 'st_marker.csv')

st_p.to_csv(save + 'st_p.csv')

scst_p.to_csv(save + 'scst_p.csv')