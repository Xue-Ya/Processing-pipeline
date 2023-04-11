import scanpy as sc
import pandas as pd

def filter(adata, Cell_type):

    invalid_groups_selected = adata.obs[Cell_type].value_counts().loc[lambda x: x < 2].index
    cell = pd.DataFrame(adata.obs[Cell_type])

    cell_filter = cell.loc[cell[Cell_type].isin(invalid_groups_selected)]

    adata = adata[~adata.obs_names.isin(cell_filter.index)]

    new_cell = cell.loc[~cell[Cell_type].isin(invalid_groups_selected)]
    adata.obs[Cell_type] = pd.Series(new_cell[Cell_type].values, index=new_cell.index)

    cell_filter = pd.Series(new_cell[Cell_type], index=new_cell.index,dtype=str)
    cell_filter = pd.Series(cell_filter, index=cell_filter.index,dtype='category')

    adata.obs[Cell_type] = cell_filter
    
    return adata