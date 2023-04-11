#!/usr/bin/env /home/xuezhengyang/.conda/envs/FYP/bin/python
import os,sys

root_dir = os.path.join(os.getcwd(),'')
if root_dir not in sys.path:
    sys.path.append(root_dir)

from src_SVI.SPIDER import *
# from src_SVI.run_somde import *
import scanpy as sc
import anndata as ad
import pandas as pd

# dataset = "Dataset1"
dataset = "Bulk"


if dataset == "Dataset1":
    regions=["st_12","st_12_2","st_14","st_14_2","st_16","st_16_2","st_17","st_17_2","st_19","st_20_1","st_20_2","st_20_3"]
    # "st_12","st_12_2","st_14",
    # ,"st_12_2","st_14","st_14_2"
    # regions=["st_17"]
    # regions=["st_16","st_16_2","st_17","st_17_2","st_19","st_20_1","st_20_2","st_20_3"]
    # "st_16","st_16_2",
    for st in regions:
        print(st)
        
        work_dict = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/st/'

        # repair_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/' + rna + '_' + st + '/Result_Repair/'
                
        output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/' + st +'/SVG/'
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        st_path = work_dict + st+'.h5ad'

        st_adata = sc.read(st_path)

        st_adata.obs['row'] = st_adata.obs['array_row']
        st_adata.obs['col'] = st_adata.obs['array_col']

        # # scST_path = repair_path + 'sprout/sc_agg_exp.tsv'
        # sc_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/sc/' + rna + '.h5ad'
        # # corr_path = repair_path + 'spot_raw_coord_best.tsv'
        # sc_cell_type = 'new_celltype'

        # sc_adata = sc.read(sc_path)

        # st_adata = st_adata[sc_adata.obs[sc_cell_type]]

        st_adata.uns['log1p']["base"] = None 
        sc.pp.highly_variable_genes(adata=st_adata, n_top_genes = 5000,subset=True)
        
        st_adata.var_names_make_unique()
        st_adata.obs_names_make_unique()
        
        spider = SPIDER()

        idata,meta_adata = SPIDER.find_svi(spider,idata = st_adata, out_f=  output_path, overwrite=False)
        idata.write_h5ad(f'{output_path}/adata.h5ad')
        meta_adata .write_h5ad(f'{output_path}/meta_adata.h5ad')

        spider.vis_pattern(idata, show_SVI=10)
        plt.tight_layout()
        plt.savefig(f'{output_path}lpattern_full.png', dpi=600,bbox_inches='tight')
        plt.close()
else:
    dorners = ['_TD1_','_TD2_','_TD3_', '_TD5_','_TD6_','_TD8_']
    # '_TD1_','_TD5_','_TD6_','_TD8_'
    # '_TD3_', '_TD5_','_TD6_','_TD8_'
    # '_TD1_','_TD2_',
    # '_TD3_',
    for dorner in dorners:
        print(dorner)
        work_dict = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/filtered/st/'

        # repair_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/' + rna + '_' + st + '/Result_Repair/'
                
        output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/bulk_dataset1/' + dorner +'/SVG/'
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        st_path = work_dict + 'st' + dorner+'.h5ad'

        st_adata = sc.read(st_path)

        st_adata.obs['row'] = st_adata.obs['array_row']
        st_adata.obs['col'] = st_adata.obs['array_col']
        
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

        new_marker = []
        marker = marker1_Epithelial + marker2_T_NK + marker3_Myeloid + marker4_B + marker5_Endothelial + marker6_Fibroblasts
        for i in st_adata.var_names:
            if i in marker:
                new_marker.append(i)
                # new_point_marker.append(i)

        new_marker = list(set(new_marker))
        
        sc.pp.log1p(st_adata)
        sc.pp.highly_variable_genes(adata=st_adata, n_top_genes = 5000)
        
        index = list(st_adata.var['highly_variable'][st_adata.var['highly_variable'] == True].index)
        
        retD = list(set(new_marker).difference(set(index)))
        
        for i in retD:
            st_adata.var['highly_variable'][i] = True
            
        st_adata = st_adata[:,st_adata.var['highly_variable']==True]
        
        st_adata.var_names_make_unique()
        st_adata.obs_names_make_unique()
        
        # st_adata.write('/home/xuezhengyang/data6/02-deconv_1/Script/notebook/test.h5ad')
        
        spider = SPIDER()

        idata,meta_adata = SPIDER.find_svi(spider,idata = st_adata, out_f=output_path, overwrite=False)
        idata.write_h5ad(f'{output_path}/adata.h5ad')
        meta_adata .write_h5ad(f'{output_path}/meta_adata.h5ad')

        spider.vis_pattern(idata, show_SVI=10)
        plt.tight_layout()
        plt.savefig(f'{output_path}lpattern_full.png', dpi=600,bbox_inches='tight')
        plt.close()