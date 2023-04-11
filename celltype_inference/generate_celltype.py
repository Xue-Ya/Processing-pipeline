import sys
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import pandas as pd
import scanpy as sc
import numpy as np
import anndata as ad
# import DeconvolutionSpot
rcParams['pdf.fonttype'] = 42
plt.style.use('default')

def filt(marker, matrix):
    # new_marker = []
    new_marker = []
    for i in matrix.var_names:
        if i in marker:
            # new_marker.append(i)
            new_marker.append(i)
            
    return new_marker

# adata_new = sc.read('/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/sc/all.h5ad')
adata_new = sc.read('/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/st/all.h5ad')


# output_path = "/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/filtered/sc"
output_path = "/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/filtered/st"


if not os.path.exists(output_path):
    os.makedirs(output_path)

# pos_filter = tissue_pos.loc[tissue_pos.iloc[:,0] == 1]

# x = pos_filter.iloc[:,3]
# y = pos_filter.iloc[:,4]

# matrix.obs['array_row'] = x
# matrix.obs['array_col'] = y

marker1_Epithelial = ["CAPS","SCGB1A1","WFDC2","KRT8","KRT19","SCGB3A1","SCGB3A2","KRT18",
          "EPCAM","MUC1","SFTPC","SFTPA1","SFTPA2","SFTPB","AGER","PGC","CAV1","CYP4B1","KRT7","ID1","SFTPD","CDH1"]

marker2_T_NK = ["GZMK","GZMB","NKG7","GNLY","KLRD1","KLRC1","TIGIT","LTB",
          "TNFRSF4","TNFRSF18","IL2RA","GZMH","CD3D","FOXP3"]

marker3_Myeloid = ["LYZ","CTSD","CCL18","CD68","CD14","HLA-DRA","HLA-DRB1","TNF",
          "FCGR3A","CD63","CD163","FCGR2A","TPSB2","TPSAB1","CPA3","HPGDS","MS4A2","AIF1","IL1RN"
          ,"CD86","THBD","S100A9","S100A8","CXCL8","FCGR3B","IL1RN","PTGS2","CD1E",
          "FCGR2A","THBD","HLA-DPB1"
          ,"CD1C","CD40"]

marker4_B = ["JCHAIN","IGHM","IGHG1","IGHA1","IGHG4","IGHA2","IGHG3","IGHG2",
          "MZB1","IGHD","CD79A","MS4A1"]

marker5_Endothelial = ["CLDN5","RAMP2","CAV1","VWF","PECAM1","CDH5","EMCN","CD34",
          "CD31","THBD"]

marker6_Fibroblasts = ["DCN","LUM","COL1A2","COL3A1","COL1A1","FN1","COL6A2","COL6A3",
          "ACTA2","COL6A1","COL5A2","COL4A2"]

point_marker = ["CAPS","SCGB1A1","GZMK","GZMB","LYZ","CTSD","JCHAIN","IGHM","CLDN5","RAMP2","DCN","LUM"]

marker = marker1_Epithelial + marker2_T_NK + marker3_Myeloid + marker4_B + marker5_Endothelial + marker6_Fibroblasts


new_marker = []
marker = marker1_Epithelial + marker2_T_NK + marker3_Myeloid + marker4_B + marker5_Endothelial + marker6_Fibroblasts
for i in adata_new.var_names:
    if i in marker:
        new_marker.append(i)
        
new_marker = list(set(new_marker))

matrix2 = adata_new[:,new_marker]

matrix_test = matrix2.copy()

sc.tl.pca(matrix_test,n_comps = 30)

#sc
# sc.pp.neighbors(matrix_test, n_neighbors=15, n_pcs=30)
# sc.tl.louvain(matrix_test,resolution=0.3,random_state=1)
#st
sc.pp.neighbors(matrix_test, n_neighbors=10, n_pcs=30)
sc.tl.louvain(matrix_test,resolution=1,random_state=1)

sc.tl.paga(matrix_test)
sc.pl.paga(matrix_test, plot=False)

#sc
# sc.tl.umap(matrix_test,random_state=1,min_dist = 0.8)

#st
sc.tl.umap(matrix_test,random_state=1)


sc.pl.umap(matrix_test,color='louvain')

sc.tl.rank_genes_groups(matrix_test, 'louvain', method='t-test')

sc.pl.rank_genes_groups(matrix_test, n_genes=5, sharey=False)

marker_gene = pd.DataFrame(matrix_test.uns['rank_genes_groups']['names']).head(5)

Epithelial = []
t_nk = []
Myeloid = []
B_cell = []
Endothelial = []
Fibroblasts = []
cell_type = [Epithelial, t_nk, Myeloid, B_cell, Endothelial, Fibroblasts ]
types = ["Epithelial", "t_nk", "Myeloid", "B_cell", "Endothelial", "Fibroblasts"]

total_marker = [filt(marker1_Epithelial,matrix_test), filt(marker2_T_NK,matrix_test), filt(marker3_Myeloid,matrix_test), filt(marker4_B,matrix_test),
                filt(marker5_Endothelial,matrix_test), filt(marker6_Fibroblasts,matrix_test)]

for marker in range(0,len(total_marker)):
    # marker_num = 0
    cell_list = cell_type[marker]
    for i in marker_gene.columns:
        count = 0
        for j in marker_gene.index:
            if marker_gene.iloc[int(j),int(i)] in total_marker[marker]:
                count = count + 1
            
            if count > 2:
                if len(cell_list) == 0:
                    cell_list.append([count, i])
                else:
                    if cell_list[-1][1] == i:
                        cell_list.pop(-1)
                        # print(cell_list[-1][1])
                    cell_list.append([count, i])
                
celltype = pd.DataFrame(index=matrix_test.to_df().index, columns=['celltype'])
for single in  cell_type:
    for cluster in single:
        index = matrix_test.obs['louvain'][matrix_test.obs['louvain'] == cluster[1]].index
        celltype.loc[index] = types[cell_type.index(single)]

path = "/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/"
matrix_test.obs['celltype'] = celltype

if celltype.isnull().sum()[0] == 0:
    
    for i in matrix_test.obs['batch'].unique():
        cell_types = pd.DataFrame(matrix_test[matrix_test.obs['batch'] == i].obs['celltype'])
        # adata = matrix_test[:,index]
        # adata.write(path +'filtered/sc/sc' + i + '.h5ad')
        cell_types.to_csv(path +'filtered/st-/st' + i + '.csv')
else:
    print(celltype.isnull().sum())
