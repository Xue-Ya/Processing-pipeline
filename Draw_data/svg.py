#!/usr/bin/env /home/xuezhengyang/.conda/envs/FYP/bin/python
import pandas as pd
import scanpy as sc
import os

#"
# regions=["st_12"]
# ,"st_12_2","st_14","st_14_2"
regions=["st_12_2","st_14","st_14_2", "st_16","st_16_2","st_17","st_17_2","st_19"]
# regions = ['st_20_3']
method = "spatalk"
dataset = "Dataset1"

deconv = ['Seurat','Tangram']

for st in regions:
    for decon in deconv:
        work_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/"+ dataset + "/" + st + "/SVG/"+decon+"/spatalk/"
        if method == "spatalk":
            out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + st + "/svg/"+decon+"/spatalk/"
        else:
            out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + st + "/svg/sprout/"
        
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        idata = sc.read(work_dir + 'meta_adata.h5ad')

        corr = idata.obs[['row','col']]

        num = idata.var.label.unique()
        num = num[num>=0]
        show_SVI = num.max() - num.min() + 1

        # svi_marker = pd.DataFrame(columns=['id'] + index)
        svi_to_plot = []
        for i in range(0,show_SVI):
            svis = idata.var[idata.var.label == i].sort_values(f'pattern_membership_{i}', ascending=False)
            svi_to_plot = svi_to_plot + svis.index.to_numpy()[:show_SVI].tolist()
            # svi_marker['pattern_membership_' + str(i)] = svi_to_plot
            # print(svi_to_plot)
            
        svi_to_plot = list(set(svi_to_plot))

        exp = idata.to_df()

        fc =  idata.obsm['pattern_score']

        index = []
        print(type(index))
        for i in range(0,show_SVI):
            index.append('pattern_membership_' + str(i))
            
        fc = pd.DataFrame(fc,index= exp.index, columns= index)

        exp = pd.concat([fc,exp],axis=1)

        exp = pd.concat([corr,exp],axis=1)


        num = idata.var.label.unique()
        num = num[num>=0]
        show_SVI = num.max() - num.min() + 1

        svi_marker = pd.DataFrame(columns= index, index=range(0,5))

        for i in range(0,5):
            svis = idata.var[idata.var.label == i].sort_values(f'pattern_membership_{i}', ascending=False)
            svi_to_plot = svis.index.to_numpy()[:show_SVI].tolist()
            
            svi_to_plot = svi_to_plot + [' '] * (5- len(svi_to_plot))
            
            svi_marker['pattern_membership_' + str(i)] = svi_to_plot
            # print(svi_to_plot)

        svi_marker.to_csv(out_dir + 'marker.csv')
        exp.to_csv(out_dir + 'exp.csv')