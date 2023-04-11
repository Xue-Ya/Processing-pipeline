import pandas as pd
import scanpy as sc
import anndata as ad



def draw_marker_genes(st,save_path,dataset,repair_path):
    data_dir = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/'+dataset+'/apart/'

    sc_cell_type = 'new_celltype'
    st_Cell_type = 'annotation'
    
    meta = pd.read_csv(repair_path+ 'sc_agg_meta.tsv',sep= '\t',index_col=1)
    
    corr = pd.read_csv(repair_path + 'spot_raw_coord_best.tsv',sep='\t',index_col=0)
    
    scst = pd.read_csv(repair_path+ 'sprout/sc_agg_exp.tsv',index_col=0)

    st = sc.read(data_dir +'st/'+ st + '.h5ad')

    scst.index = scst['index']
    scst = scst.drop(['index'],axis=1)

    # celltype = set(sc_adata.obs[sc_cell_type].index)

    cell = scst.index
    index = pd.DataFrame(index=scst.index)

    scst_adata = ad.AnnData(scst,obs = index, dtype='int32')

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

    scst_marker_gene = pd.DataFrame(scst_adata_new.uns['rank_genes_groups']['names']).head(5)
    
    #combine gene expression
    exp = scst_adata_new.to_df()
    
    exp = pd.concat([meta.loc[:,['x','y']],exp], axis=1)
    exp.loc[:,'x'] = corr[:,0]
    exp.y = corr[:,1]
    
    genes = ['x','y']
    for row in range(0,len(scst_marker_gene.index)):
        for col in range(0,len(scst_marker_gene.columns)):
            gene = scst_marker_gene.iloc[row,col]
            genes.append(gene)

    genes2 = list(set(genes))

    genes2.sort(key=genes.index)
    scst_exp = exp.loc[:,genes2]   

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
    
    #combine gene expression
    pois_key = ['array_row','array_col']
    
    exp = st_new.to_df()
    
    row = pd.DataFrame(st_new.obs[pois_key[0]])
    col = pd.DataFrame(st_new.obs[pois_key[1]])

    corr = pd.concat([row,col], axis=1)
    
    exp = pd.concat([corr,exp], axis=1)
    
    genes = ['array_row','array_col']
    for row in range(0,len(st_marker_gene.index)):
        for col in range(0,len(st_marker_gene.columns)):
            gene = st_marker_gene.iloc[row,col]
            genes.append(gene)

    
    genes2 = list(set(genes))

    genes2.sort(key=genes.index)
    st_exp = exp.loc[:,genes2]    

    result = st_new.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    st_p = pd.DataFrame(
        {group + ' '+ key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(5)

    ########## To csv
    
    import os
    
    if not os.path.exists(save_path + 'scst/'):
        os.makedirs(save_path + 'scst/')
        
    if not os.path.exists(save_path + 'st/'):
        os.makedirs(save_path + 'st/')
    
    na = scst_exp[scst_exp.isnull().T.any()].index
    scst_exp = scst_exp.drop(na)
    
    na = scst_exp[st_exp.isnull().T.any()].index
    scst_exp = st_exp.drop(na)
    
    scst_marker_gene.index.name='id'
    st_marker_gene.index.name='id'
    st_p.index.name='id'
    scst_p.index.name='id'
    scst_exp.index.name='id'
    st_exp.index.name='id'
    
    scst_marker_gene.to_csv(save_path + 'scst/scst_marker.csv')

    st_marker_gene.to_csv(save_path + 'st/st_marker.csv')

    st_p.to_csv(save_path + 'st/st_p.csv')

    scst_p.to_csv(save_path + 'scst/scst_p.csv')
    
    scst_exp.to_csv(save_path + 'scst/scst_exp.csv')
    
    st_exp.to_csv(save_path + 'st/st_exp.csv')