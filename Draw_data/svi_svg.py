import scanpy as sc
sc.settings.verbosity = 3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
# %load_ext autoreload
# %reload_ext autoreload
# %autoreload 2
root_dir = os.path.join(os.getcwd(), '..')
if root_dir not in sys.path:
    sys.path.append(root_dir)
import json
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
import anndata 
from scipy import stats
from scipy import spatial
import umap, umap.plot
import gseapy
import seaborn as sns

def draw_diff(adata, idata, title, is_human=True):
    from matplotlib_venn import venn2,venn2_circles

    organism = 'Human' if is_human else 'Mouse'
    pathway_db = 'KEGG_2021_Human' if is_human else 'KEGG_2019_Mouse'

    print(organism, pathway_db)

    # genelist = adata.uns['moranI'][:50].index.to_list()
    genelist = adata.uns['moranI'][:200].index.to_list()
    enr_res_gene = gseapy.enrichr(gene_list=genelist,
                            cutoff=0.01,
                            organism=organism,
                            gene_sets=pathway_db)
    gene_pw_list = enr_res_gene.res2d.Term.to_numpy()

    interaction_list = idata.uns['moranI'][:200].index.to_list()
    gene_list_lri = np.unique(np.concatenate([x.split('_') for x in interaction_list])).tolist()
    enr_res = gseapy.enrichr(gene_list=gene_list_lri,
                                cutoff=0.01,
                                organism=organism,
                                gene_sets=pathway_db)
    lri_pw_list = enr_res.res2d.Term.to_numpy()

    fig = plt.figure(figsize=(10, 5))
    # fig.suptitle(title, fontsize=15)
    plt.subplot(1, 3, 1)
    # moran_df = adata.uns['moranI'][:50]
    moran_df = adata.uns['moranI'][:200]
    moran_df['label'] = 'SVG'
    moran_df_svi = adata.uns['moranI'].loc[np.intersect1d(gene_list_lri,adata.uns['moranI'].index)]
    moran_df_svi['label'] = 'SVI genes'
    df = pd.concat([moran_df, moran_df_svi])
    print(moran_df.I.mean(), moran_df_svi.I.mean())
    sns.boxplot(data=df, x="label", y="I", order=["SVI genes", "SVG"], palette={"SVI genes":"#098154", "SVG": "#c72e29"})
    plt.title('Moran I')
    
    plt.subplot(1, 3, 2)
    g=venn2(subsets = [set(gene_list_lri),set(genelist)], #绘图数据集
        set_labels = ('SVI genes', 'SVG'), #设置组名
        set_colors=("#098154","#c72e29"),#设置圈的颜色，中间颜色不能修改
        alpha=0.6,#透明度
        normalize_to=1.0,#venn图占据figure的比例，1.0为占满
       )
    plt.title('Gene Overlap')
    plt.subplot(1, 3, 3)
    g=venn2(subsets = [set(lri_pw_list),set(gene_pw_list)], #绘图数据集
        set_labels = ('SVI genes', 'SVG'), #设置组名
        set_colors=("#098154","#c72e29"),#设置圈的颜色，中间颜色不能修改
        alpha=0.6,#透明度
        normalize_to=1.0,#venn图占据figure的比例，1.0为占满
       )
    plt.title('Enriched Pathway Overlap')

regions=["st_12","st_12_2","st_14","st_14_2", "st_16","st_16_2","st_17","st_17_2","st_19"]
# regions = ['st_12']
method = "spatalk"
dataset = "Dataset1"

deconv = ['Seurat','Tangram']

for st in regions:
    for decon in deconv:
        if method == "spatalk":
            idata = sc.read("/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/SVI/scST/"+decon +'_'+ st + "/meta_idata.h5ad")
            adata = sc.read("/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/"+ dataset + "/" + st + "/SVG/"+decon+"/spatalk/meta_adata.h5ad")
            out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + st + "/svg-svi/"+decon+"/spatalk/"
        else:
            out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + st + "/svg/sprout/"
            
        
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            
        is_human=True
        from matplotlib_venn import venn2,venn2_circles

        organism = 'Human' if is_human else 'Mouse'
        pathway_db = 'KEGG_2021_Human' if is_human else 'KEGG_2019_Mouse'

        print(organism, pathway_db)

        # genelist = adata.uns['moranI'][:50].index.to_list()
        genelist = adata.uns['moranI'][:200].index.to_list()
        enr_res_gene = gseapy.enrichr(gene_list=genelist,
                                cutoff=0.01,
                                organism=organism,
                                gene_sets=pathway_db)
        gene_pw_list = enr_res_gene.res2d.Term.to_numpy()

        interaction_list = idata.uns['moranI'][:200].index.to_list()
        gene_list_lri = np.unique(np.concatenate([x.split('_') for x in interaction_list])).tolist()
        enr_res = gseapy.enrichr(gene_list=gene_list_lri,
                                    cutoff=0.01,
                                    organism=organism,
                                    gene_sets=pathway_db)
        lri_pw_list = enr_res.res2d.Term.to_numpy()

        fig = plt.figure(figsize=(10, 5))
        # fig.suptitle(title, fontsize=15)
        plt.subplot(1, 3, 1)
        # moran_df = adata.uns['moranI'][:50]
        moran_df = adata.uns['moranI'][:200]
        moran_df['label'] = 'SVG'
        moran_df_svi = adata.uns['moranI'].loc[np.intersect1d(gene_list_lri,adata.uns['moranI'].index)]
        moran_df_svi['label'] = 'SVI genes'
        df = pd.concat([moran_df, moran_df_svi])
        print(moran_df.I.mean(), moran_df_svi.I.mean())
        sns.boxplot(data=df, x="label", y="I", order=["SVI genes", "SVG"], palette={"SVI genes":"#098154", "SVG": "#c72e29"})
        plt.title('Moran I')

        plt.subplot(1, 3, 2)
        g=venn2(subsets = [set(gene_list_lri),set(genelist)], #绘图数据集
            set_labels = ('SVI genes', 'SVG'), #设置组名
            set_colors=("#098154","#c72e29"),#设置圈的颜色，中间颜色不能修改
            alpha=0.6,#透明度
            normalize_to=1.0,#venn图占据figure的比例，1.0为占满
            )
        plt.title('Gene Overlap')
        plt.subplot(1, 3, 3)
        g=venn2(subsets = [set(lri_pw_list),set(gene_pw_list)], #绘图数据集
            set_labels = ('SVI genes', 'SVG'), #设置组名
            set_colors=("#098154","#c72e29"),#设置圈的颜色，中间颜色不能修改
            alpha=0.6,#透明度
            normalize_to=1.0,#venn图占据figure的比例，1.0为占满
            )
        plt.title('Enriched Pathway Overlap')
        
        
        #  box plot
        box = df.drop(['pval_norm','var_norm','pval_z_sim','pval_sim','var_sim','pval_norm_fdr_bh','pval_z_sim_fdr_bh','pval_sim_fdr_bh'],axis=1)

        box.index.name='id'

            
        box.to_csv(out_dir + 'box_plot.csv')
        # lri_pw_list
        
        #Venn 1
        venn_left = pd.DataFrame(index = list(set(gene_list_lri)),columns=['label'])
        venn_left.loc[:,'label'] = 'SVI'
        svg = pd.DataFrame(index = list(set(genelist)),columns=['label'])
        svg.loc[:,'label'] = 'SVG'

        venn_left = venn_left.append(svg)

        venn_left.index.name='id'

        retD = list(set(gene_list_lri).intersection(set(genelist)))
        venn_left.loc[retD] = 'Gene Overlap'

        venn_left = venn_left[~venn_left.index.duplicated(keep='first')]

        venn_left.to_csv(out_dir +'gene_overlap.csv')
        # gene_list_lri = set(gene_list_lri)
        # venn_left[gene_list_lri]['label'] = 'SVI'
        # svi = gene_list_lri
        
        
        # Venn 2
        venn_right = pd.DataFrame(index = list(set(gene_pw_list)),columns=['label'])
        venn_right.loc[:,'label'] = 'SVI'
        svg = pd.DataFrame(index = list(set(lri_pw_list)),columns=['label'])
        svg.loc[:,'label'] = 'SVG'

        venn_right = venn_right.append(svg)

        retD = list(set(gene_pw_list).intersection(set(lri_pw_list)))
        venn_right.loc[retD] = 'Enriched Pathway Overlap'

        venn_right.index.name='id'

        venn_right = venn_right[~venn_right.index.duplicated(keep='first')]

        venn_right.to_csv(out_dir +'pathway_overlap.csv')

        # lri_pw_list
        
        moran_df_svi = adata.uns['moranI'].loc[np.intersect1d(gene_list_lri,adata.uns['moranI'].index)]
        moran_df_svi['label'] = 'SVI genes'
        # df = pd.concat([moran_df, moran_df_svi])
        # moran_df_svi
        
        res = []
        res.append([x.split('_') for x in interaction_list])
        res = res[0]
        
        row = []
        row.append([x[0] for x in res])
        row = row[0]

        col = []
        col.append([x[1] for x in res])
        col = col[0]
        
        # Chord
        chord = pd.DataFrame(index=range(0,len(col)),columns=['Col','Value','PValue','Metric','Method','Feature','Group'])
        
        expression = idata.to_df()
        exp = expression.loc[:,interaction_list].mean()
        
        chord.index = row
        chord['Col'] = col
        chord['Value'] = exp.values
        # ,quanTIseq,c_tumor_stage,Stage I
        
        chord['PValue'] = 0
        chord['Metric'] = 'Pearson Correlation'
        chord['Method'] = 'quanTIseq'
        chord['Feature'] = 'c_tumor_stage'
        chord['Group'] = 'Stage I'
        
        chord = chord[~chord.index.duplicated(keep='first')]
        chord.index.name='ROW'

        chord.to_csv(out_dir + 'chord.csv')