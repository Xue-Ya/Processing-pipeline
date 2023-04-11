import pandas as pd
import scanpy as sc
import anndata as ad

def spatalk(st,deconv,dataset,save_path):
    
    dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/'+dataset+'/' + st + '/Result_Repair/'+deconv+'/spatalk/'
    # dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/stage2_st_12_2/Result_Repair/spatalk/'
    exp_path = dir + 'exp.csv'
    meta_path = dir + 'meta.csv'
    genes_path = dir + 'genes.csv'
    
    meta =pd.read_csv(meta_path,index_col=0)
    genes =pd.read_csv(genes_path)
    exp = pd.read_csv(exp_path)

    for i in range(0,meta['celltype'].size):
        char = meta['celltype'].iloc[i]
        char = char.replace('_', ' ')
        meta['celltype'].iloc[i] = char
        
    #Expression
    exp.index = genes.values[:,0]
    exp = exp.T
    
    celltype  = meta.celltype

    # for i in range(0,celltype.index.size):
        # new_cell = celltype.iloc[i][1:]
        # celltype.iloc[i] = new_cell
    
    #Meta
    celltype = meta.loc[:,['y','x','celltype']]
      
    
    scst_adata = ad.AnnData(exp,obs = meta, dtype='int32')
    
    corr = meta.loc[:,['x','y']]
    
    scst_adata.obsm['spatial'] = corr
    
    Cell_type = 'celltype'
    
    # Drop samples only show once
    invalid_groups_selected = scst_adata.obs[Cell_type].value_counts().loc[lambda x: x < 2].index
    cell = pd.DataFrame(scst_adata.obs[Cell_type])

    cell_filter = cell.loc[cell[Cell_type].isin(invalid_groups_selected)]

    scst_adata = scst_adata[~scst_adata.obs_names.isin(cell_filter.index)]

    new_cell = cell.loc[~cell[Cell_type].isin(invalid_groups_selected)]
    scst_adata.obs[Cell_type] = pd.Series(new_cell[Cell_type].values, index=new_cell.index)

    cell_filter = pd.Series(new_cell[Cell_type], index=new_cell.index,dtype=str)
    cell_filter = pd.Series(cell_filter, index=cell_filter.index,dtype='category')

    scst_adata.obs[Cell_type] = cell_filter
    
    
    #rank genes
    
    sc.pp.log1p(scst_adata)

    sc.tl.rank_genes_groups(scst_adata, 'celltype',  method='wilcoxon')

    sc.pl.rank_genes_groups(scst_adata,sharey=False,n_genes=10,show = False)
    
    scst_marker_gene = pd.DataFrame(scst_adata.uns['rank_genes_groups']['names']) .head(5)
    
    #combine gene expression
    exp = scst_adata.to_df()
    
    exp = pd.concat([meta.loc[:,['x','y']],exp], axis=1)
    # exp.loc[:,'x'] = corr[:,0]
    # exp.y = corr[:,1]
    
    genes = ['x','y']
    for row in range(0,len(scst_marker_gene.index)):
        for col in range(0,len(scst_marker_gene.columns)):
            gene = scst_marker_gene.iloc[row,col]
            genes.append(gene)

    genes2 = list(set(genes))

    genes2.sort(key=genes.index)
    scst_exp = exp.loc[:,genes2]   

    result = scst_adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    scst_p = pd.DataFrame(
        {group + ' '+ key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(5)
    
    
    scst_p.index.name='id'
    scst_exp.index.name='id'
    celltype.index.name='id'
    scst_marker_gene.index.name='id'
     
    na = scst_exp[scst_exp.isnull().T.any()].index
    scst_exp = scst_exp.drop(na)
    
    scst_marker_gene.to_csv(save_path + '/marker_gene/'+ deconv + '/spatalk/scst_marker.csv')
    
    scst_p.to_csv(save_path + '/marker_gene/'+ deconv + '/spatalk/scst_p.csv')
    
    scst_exp.to_csv(save_path + '/marker_gene/'+ deconv + '/spatalk/scst_exp.csv')
    
    #Draw scst spot
    celltype.to_csv(save_path + '/scst/'+ deconv + '/scst_spot_spatalk.csv')