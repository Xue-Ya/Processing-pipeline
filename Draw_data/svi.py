#!/usr/bin/env /home/xuezhengyang/.conda/envs/FYP/bin/python
import pandas as pd
import scanpy as sc
import os



# 
# regions=["st_12_2","st_14","st_14_2"]
# "st_12","st_12_2","st_14","st_14_2"

# error
# Tangram: "st_16","st_16_2", "st_17_2"

regions=["st_17","st_17_2"]
# regions = ['st_20_3']
method = "spatalk"
dataset = "dataset1"
deconv = ['Seurat_','Tangram_']

for st in regions:
    for decon in  deconv:
    # /home/xuezhengyang/data6/02-deconv_1/Script/FigureData/"+dataset+"/stage2_st_12_2/SVI/sprout/stage2_st_12_2_sprout/idata.h5ad
        if method == "spatalk":
            work_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/SVI/scST/"+decon + st + "/"
            out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + st + "/svi/"+decon+"/spatalk/"
        else:
            work_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/"+dataset+"/" + st + "/SVI/sprout/"
            out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + st + "/svi/sprout/"
        
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        idata = sc.read(work_dir + 'meta_idata.h5ad')

        corr = idata.obs[['row','col']]

        num = idata.var.label.unique()
        num = num[num>=0]
        show_SVI = num.max() - num.min() + 1
        # show_SVI = idata.var.label.max() - idata.var.label.min() 

        print(show_SVI)
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
            print(index)
            
        fc = pd.DataFrame(fc,index= exp.index, columns= index)

        exp = pd.concat([fc,exp],axis=1)

        exp = pd.concat([corr,exp],axis=1)


        show_SVI = show_SVI

        svi_marker = pd.DataFrame(columns= index, index=range(0,show_SVI))

        for i in range(0,show_SVI):
            svis = idata.var[idata.var.label == i].sort_values(f'pattern_membership_{i}', ascending=False)
            svi_to_plot = svis.index.to_numpy()[:show_SVI].tolist()
            
            svi_to_plot = svi_to_plot + [' '] * (show_SVI- len(svi_to_plot))
            
            svi_marker['pattern_membership_' + str(i)] = svi_to_plot
            # print(svi_to_plot)

        svi_marker.to_csv(out_dir + 'marker.csv')
        exp.to_csv(out_dir + 'exp.csv')