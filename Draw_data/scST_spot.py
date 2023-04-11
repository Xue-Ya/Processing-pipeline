import pandas as pd

def draw_scst_spot(st,save_path,repair_path):

    meta = pd.read_csv(repair_path+ 'sc_agg_meta.tsv',sep= '\t',index_col=0)
    corr = pd.read_csv(repair_path + 'spot_raw_coord_best.tsv',sep='\t',index_col=0)

    final = meta.loc[:,['x','y','celltype']]
    final.loc[:,'x'] = corr.iloc[:,1]
    final.y = corr.iloc[:,0]
    final.index = meta.loc[:,'spot']
    
    final.index.name='id'
    final.to_csv(save_path + '/scst_spot_sprout.csv')