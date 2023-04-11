import pandas as pd
import scanpy as sc
import anndata as ad

regions = ["st_12","st_12_2","st_14","st_14_2","st_16","st_16_2","st_17","st_17_2"]
dataset = "Dataset1"

for st in regions:
    st_Cell_type = 'annotation'

    data_dir = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/'+dataset+'/apart/'

    sc_cell_type = 'new_celltype'
    st_Cell_type = 'annotation'

    save_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Figure/' + st + '/'

    st = sc.read(data_dir +'st/'+ st + '.h5ad')

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

    st_marker_gene.dropna()
    st_p.dropna()
    st_exp.dropna()

    st_marker_gene.index.name='id'
    st_p.index.name='id'
    st_exp.index.name='id'
    
    st_exp.index = st_exp.index.str.replace('-','_')

    st_marker_gene.to_csv(save_path + 'st/st_marker.csv')

    st_p.to_csv(save_path + 'st/st_p.csv')

    st_exp.to_csv(save_path + 'st/st_exp.csv')