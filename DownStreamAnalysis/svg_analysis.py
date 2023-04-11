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

# rna = sys.argv[1]
# st = sys.argv[2]

# rna="stage2"
# rna="stage3"
# regions=["st_12","st_12_2","st_14","st_14_2"]
# "st_12","st_12_2","st_14",
# ,"st_12_2","st_14","st_14_2"
# regions=["st_17"]
# regions = ["st_12","st_12_2","st_14","st_14_2" "st_16","st_16_2","st_17","st_16","st_16_2","st_17_2","st_19"]
st = sys.argv[1]
dorner = sys.argv[1]
# ,"st_12_2","st_14","st_14_2" "st_16","st_16_2","st_17","st_16","st_16_2","st_17_2","st_19","st_20_1","st_20_2","st_20_3"
# "st_16","st_16_2",

# dataset = "Dataset1"
dataset = "bulk_dataset1"
# dataset = "bulk"
# method = "sprout"
deconv = ['Tangram', 'Seurat']
method = "spatalk"

if dataset == "Dataset1":
    for decon in deconv:
        if method == "spatalk":
            # for st in regions:
            print("info!:", st, dataset, decon)
            
            work_dict = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/'  + st + '/Result_Repair/'+decon+'/spatalk/'

            repair_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/' + st + '/Result_Repair/'
                    
            output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/' + st +'/SVG/'+ decon + '/spatalk/'
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            scST_path = work_dict+st+'.h5ad'

            scst_adata = sc.read(scST_path)

            scst_adata.obs['row'] = scst_adata.obs['x']
            scst_adata.obs['col'] = scst_adata.obs['y']
            
            sc.pp.highly_variable_genes(adata=scst_adata, n_top_genes = 5000,subset=True)
            
            # print(scst_adata.shape)
            scst_adata.var_names_make_unique()
            scst_adata.obs_names_make_unique()
            
            spider = SPIDER()

            idata,meta_adata = SPIDER.find_svi(spider,idata = scst_adata, out_f=  output_path, overwrite=False)
            idata.write_h5ad(f'{output_path}/adata.h5ad')
            meta_adata .write_h5ad(f'{output_path}/meta_adata.h5ad')

            spider.vis_pattern(idata, show_SVI=10)
            plt.tight_layout()
            plt.savefig(f'{output_path}lpattern_full.png', dpi=600,bbox_inches='tight')
            plt.close()
        else:
            print(st, dataset, method)
            work_dict = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/' + st + '/Result_Repair/sprout/'
                    
            output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/' + st +'/SVG/sprout/'
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            scST_path = work_dict +st+'_sprout.h5ad'

            scst_adata = sc.read(scST_path)

            scst_adata.obs['row'] = scst_adata.obs['x']
            scst_adata.obs['col'] = scst_adata.obs['y']

            # scST_path = repair_path + 'sprout/sc_agg_exp.tsv'
            # corr_path = repair_path + 'spot_raw_coord_best.tsv'
            # scst_adata = scst_adata[sc_adata.obs[sc_cell_type]]

            spider = SPIDER()

            idata,meta_adata = SPIDER.find_svi(spider,idata = scst_adata, out_f=  output_path, overwrite=False,threshold=0.3)
            idata.write_h5ad(f'{output_path}/adata.h5ad')
            meta_adata .write_h5ad(f'{output_path}/meta_adata.h5ad')

            spider.vis_pattern(idata, show_SVI=10)
            plt.tight_layout()
            plt.savefig(f'{output_path}lpattern_full.png', dpi=600,bbox_inches='tight')
            plt.close()
else:
    #,'_TD5_' ,'_TD6_','_TD8_'
    if method == "spatalk":
        
        for decon in deconv:
            print(dorner, dataset, method)
            work_dict = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/bulk_dataset1/' + dorner + '/Result_Repair/'+decon+'/'
                    
            output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/bulk_dataset1/'+dorner +'/SVG/' + decon + '/spatalk/'
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            scST_path = work_dict + 'spatalk/'+dorner+'.h5ad'

            scst_adata = sc.read(scST_path)

            scst_adata.obs['row'] = scst_adata.obs['x']
            scst_adata.obs['col'] = scst_adata.obs['y']

            # # scST_path = repair_path + 'sprout/sc_agg_exp.tsv'
            # sc_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/sc/' + rna + '.h5ad'
            # # corr_path = repair_path + 'spot_raw_coord_best.tsv'
            # sc_cell_type = 'new_celltype'

            # sc_adata = sc.read(sc_path)

            # scst_adata = scst_adata[sc_adata.obs[sc_cell_type]]
            # print(scst_adata.shape)
            
            sc.pp.highly_variable_genes(adata=scst_adata, n_top_genes = 5000,subset=True)
            
            # print(scst_adata.shape)
            
            spider = SPIDER()

            idata,meta_adata = SPIDER.find_svi(spider,idata = scst_adata, out_f=  output_path, overwrite=True)
            idata.write_h5ad(f'{output_path}/adata.h5ad')
            meta_adata .write_h5ad(f'{output_path}/meta_adata.h5ad')

            spider.vis_pattern(idata, show_SVI=10)
            plt.tight_layout()
            plt.savefig(f'{output_path}lpattern_full.png', dpi=600,bbox_inches='tight')
            plt.close()
    else:
        print(dorner, dataset, method)
        work_dict = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/bulk_dataset1/' + dorner + '/Result_Repair/sprout/'
                
        output_path = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/bulk_dataset1/'+ dorner +'/SVG/sprout/'
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        scST_path = work_dict +dorner+'_sprout.h5ad'

        scst_adata = sc.read(scST_path)

        scst_adata.obs['row'] = scst_adata.obs['x']
        scst_adata.obs['col'] = scst_adata.obs['y']

        # scST_path = repair_path + 'sprout/sc_agg_exp.tsv'
        # corr_path = repair_path + 'spot_raw_coord_best.tsv'
        # scst_adata = scst_adata[sc_adata.obs[sc_cell_type]]

        spider = SPIDER()

        idata,meta_adata = SPIDER.find_svi(spider,idata = scst_adata, out_f=  output_path, overwrite=False,threshold=0.3)
        idata.write_h5ad(f'{output_path}/adata.h5ad')
        meta_adata .write_h5ad(f'{output_path}/meta_adata.h5ad')

        spider.vis_pattern(idata, show_SVI=10)
        plt.tight_layout()
        plt.savefig(f'{output_path}lpattern_full.png', dpi=600,bbox_inches='tight')
        plt.close()