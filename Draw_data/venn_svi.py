#!/usr/bin/env /home/xuezhengyang/.conda/envs/FYP/bin/python
import pandas as pd
import scanpy as sc
import os
import numpy as np

threshold = 0.01
dataset = "Dataset1"
# dataset = "bulk_dataset1"
method = "spatalk"
# method = ""

import matplotlib.pyplot as plt
import seaborn as sns


i_or_g = "SVI"

deconv = ['Seurat_','Tangram_']

def vis_pattern(idata, label='', pos_key=[], obsm_key='pattern_score', traj_coords=[], show_SVI=0):
    plt.figure(figsize=(30, 50))
    if label != '':
        plt.subplot(20, 11, 1)
        ax=sns.scatterplot(data=idata.uns['cell_meta'], x=pos_key[0], y=pos_key[1], hue=label, s=10, palette ='tab20')
        # ax=sns.scatterplot(data=idata.obs, x='row', y='col', hue=label, s=10, palette ='tab20')
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1.05),ncol=int(len(idata.uns['cell_meta'][label].unique())/8+1), title=None, frameon=False)
        if len(traj_coords) != 0:
            for traj_cood in traj_coords:
                plt.quiver(traj_cood[0],traj_cood[2],traj_cood[1]-traj_cood[0],traj_cood[3]-traj_cood[2],color = 'red', angles='xy', scale_units='xy', scale = 1)
        plt.axis('equal')
        plt.axis('off')
        base = 9
    else:
        base = 1
    if show_SVI==0:
        for i in range(idata.obsm[obsm_key].shape[1]):
            plt.subplot(20, 11, i + base)
            im=plt.scatter(idata.obs['row'], idata.obs['col'], c=idata.obsm[obsm_key][:,i], s=1, cmap='plasma')
            if len(traj_coords) != 0:
                for traj_cood in traj_coords:
                    plt.quiver(traj_cood[0],traj_cood[2],traj_cood[1]-traj_cood[0],traj_cood[3]-traj_cood[2],color = 'red', angles='xy', scale_units='xy', scale = 1)
            plt.axis('equal')
            plt.axis('off')
            plt.title('Pattern {} - {} LRIs'.format(i, np.sum(idata.var.label == i)) )
            plt.colorbar(im,fraction=0.046, pad=0.04)
    else:
        for i in range(idata.obsm[obsm_key].shape[1]):
            plt.subplot(20, 11, base)
            im=plt.scatter(idata.obs['row'], idata.obs['col'], c=idata.obsm[obsm_key][:,i], s=1, cmap='plasma')
            if len(traj_coords) != 0:
                for traj_cood in traj_coords:
                    plt.quiver(traj_cood[0],traj_cood[2],traj_cood[1]-traj_cood[0],traj_cood[3]-traj_cood[2],color = 'red', angles='xy', scale_units='xy', scale = 1)
            plt.axis('equal')
            plt.axis('off')
            plt.title('Pattern {} - {} LRIs'.format(i, np.sum(idata.var.label == i)))
            plt.colorbar(im,fraction=0.046, pad=0.04)
            base += 1

            svis = idata.var[idata.var.label == i].sort_values(f'pattern_membership_{i}', ascending=False)
            svi_to_plot = svis.index.to_numpy()[:show_SVI]
            memberships = svis[f'pattern_membership_{i}'].to_numpy()[:show_SVI]
            for j in range(show_SVI):
                if j < len(svi_to_plot):
                    plt.subplot(20,11, base)
                    im=plt.scatter(idata.obs['row'], idata.obs['col'], c=idata.to_df()[svi_to_plot[j]], s=1, cmap='plasma')
                    if len(traj_coords) != 0:
                        for traj_cood in traj_coords:
                            plt.quiver(traj_cood[0],traj_cood[2],traj_cood[1]-traj_cood[0],traj_cood[3]-traj_cood[2],color = 'red', angles='xy', scale_units='xy', scale = 1)
                    plt.axis('equal')
                    plt.axis('off')
                    plt.colorbar(im,fraction=0.046, pad=0.04)
                    plt.title(f'{svi_to_plot[j]}')
                base += 1



if dataset == "Dataset1":
    # regions=["st_12_2"]
    # ,"st_12_2","st_14","st_14_2"
    # ,"st_20_1","st_20_2","st_20_3"
    regions=["st_12_2","st_14","st_14_2","st_16","st_16_2","st_17","st_17_2"]
    # regions = ['st_20_3']
    for st in regions:
        
        for decon in deconv:
            print(st)
            if method == "spatalk":
                work_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/SVI/scST/"+decon + st + "/"
                # dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/' + st + '/Result_Repair/spatalk/'
                out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + st + "/venn/Tangram_spatalk/"+i_or_g+"/"
            elif method == "sprout":
                work_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/SVI/scST/"+decon + st + "/"
                dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/' + st + '/Result_Repair/sprout/'
                out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + st + "/venn/Tangram_sprout/"+i_or_g+"/"
            
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            idata = sc.read(work_dir + '/meta_idata.h5ad')
            # idata = sc.read("/data5/lishiying/lung_data/stage2_st_12/meta_idata.h5ad")
            
            label = idata.var['label'].copy()
            
            label[label>-1] = 1
            
            # idata.var['is_svi'] = label.copy()
            
            idata.var[[f'pattern_correlation_{x}' for x in range(idata.obsm['pattern_score'].shape[1])]] = 0
            corr_df=pd.concat([idata[:,idata.var['is_svi']==1].to_df(),pd.DataFrame(idata.obsm['pattern_score'],index=idata.obs_names)],axis=1).corr().loc[idata[:,idata.var['is_svi']==1].var_names, range(idata.obsm['pattern_score'].shape[1])]
            idata.var.loc[idata[:,idata.var['is_svi']==1].var_names, [f'pattern_correlation_{x}' for x in range(idata.obsm['pattern_score'].shape[1])]] = corr_df.to_numpy()
            # svis = idata.var[idata.var.label == i].sort_values(f'pattern_correlation_{i}', ascending=False)
            
            df_corr = pd.DataFrame()
            
            for i in range(0,idata.obsm['pattern_score'].shape[1]):
                df_corr = pd.concat([df_corr,idata.var[f'pattern_correlation_{i}']],axis=1)
            
            exp = idata.to_df()
            
            # if method == "spatalk" or method == "sprout":
            #     exp_path = dir + 'exp.csv'
            #     genes_path = dir + 'genes.csv'
                
            #     exp = pd.read_csv(exp_path)
            #     genes =pd.read_csv(genes_path)
                
            #     exp.index = genes.values[:,0]
            #     exp = exp.T
            # else:
            #     exp = sc.read(dir)
            #     exp.obs_names_make_unique()
            #     exp.var_names_make_unique()
            #     exp = exp.to_df()
            # idata.uns['moranI'] = pd.read_csv(work_dir+'moranI.csv',index_col=0)
            # idata.uns['gearyC'] = pd.read_csv(work_dir+'gearyC.csv',index_col=0)
            # idata.uns['nnSVG'] = pd.read_csv(work_dir+'nnSVG.csv',index_col=0)
            
            
            I_value = idata.uns['moranI']["I"].dropna()
            
            C_value = idata.uns['gearyC']["C"].dropna()
            
            filt = list(I_value.index)
            
            
            idata_true = idata[:,filt]
            # idata_true
            
            # exp = exp.loc[idata.obs_names,idata.var_names]
            
            index = idata_true.var_names
            
            SOMDE_index = idata_true.uns['SOMDE'].index
            SOMDE_index = index[SOMDE_index]
            SOMDE_index = list(SOMDE_index)
            
            SPARKX_df = idata_true.uns['SPARKX']['adjustedPval']
            SPARKX_df = SPARKX_df[(SPARKX_df<threshold)]
            SPARKX_index = list(SPARKX_df.index)
            if len(SPARKX_index) == 0:
                SPARKX_df = idata_true.uns['SPARKX']['combinedPval']
                SPARKX_df = SPARKX_df[(SPARKX_df<threshold)]
                SPARKX_index = list(SPARKX_df.index)
                
            SpatialDE_df = idata_true.uns['SpatialDE']['padj']
            SpatialDE_df = SpatialDE_df[(SpatialDE_df<threshold)]
            SpatialDE_index = list(SpatialDE_df.index)
            
            nnSVG_index = idata_true.uns['nnSVG'].index
            # nnSVG_index = index[nnSVG_index]
            nnSVG_index = list(nnSVG_index)
            
            if 'scGCO' in idata.uns_keys():
                idata.uns['scGCO'].index = idata.uns['scGCO']["Unnamed: 0"]
                scGCO_df = idata.uns['scGCO']['p_value']
                values = idata.uns['scGCO']['p_value'].values
                for i in range(0,len(values)):
                    values[i] = values[i][1:-1]
                    values[i] = values[i].split(',')
                    values[i] = float(values[i][0])
                idata.uns['scGCO']['p_value'] = values
                scGCO_df = scGCO_df[(scGCO_df<threshold)]
                scGCO_index = list(scGCO_df.index)
                
                if len(scGCO_index) == 0:
                    scGCO_df = idata.uns['nnSVG']['combinedPval']
                    scGCO_df = scGCO_df[(scGCO_df<threshold)]
                    scGCO_index = list(scGCO_df.index)
            
            else:
                scGCO_index = []
        
            
            
            df = pd.DataFrame(index=index, columns=["Label","SOMDE", "SPARKX", "SpatialDE","nnSVG","scGCO", "I", "C","Corr"])
            
            df[["Label","SOMDE", "SPARKX", "SpatialDE","nnSVG","scGCO", "I", "C","Corr"]] = 0
            
            df.loc[SOMDE_index,"SOMDE"] = 1
            
            df.loc[nnSVG_index,"nnSVG"] = 1
            
            df.loc[SPARKX_index,"SPARKX"] = 1
            
            # df.loc[scGCO_index,"scGCO"] = 1
            
            df.loc[SpatialDE_index,"SpatialDE"] = 1
            
            df["Label"] = label
            
            # print(df_corr)
            
            # df_corr.dropna(how='any',axis=0)
            
            # df_corr.to_csv('df.csv')
            
            df["Corr"] = df_corr.idxmax(axis=1)
            
            # df["Corr"].dropna(how='any',axis=0)
            
            # df["Corr"].to_csv('corr.csv')
            
            for i in range(0,df["Corr"].size):
                
                # print(df["Corr"][i])
                sp = df["Corr"][i].split('_')
                sp[1] = 'membership'
                res = sp[0] + '_' + sp[1] + '_' + sp[2]
                df["Corr"][i] = res
                    
            
            df["C"] = C_value
            df["I"] = I_value
            
            corr = idata.obs[['row','col']]
            
            exp = pd.concat([corr,exp],axis=1)
            
            exp.index.name='id'
            df.index.name='id'
            
            df = df.sort_values(by="I" , ascending=False)
            
            exp.to_csv(out_dir + "exp.csv")
            df.to_csv(out_dir + "meta.csv")
            
            vis_pattern(idata, show_SVI=10)
            plt.tight_layout()
            plt.savefig(f'{work_dir}lpattern_full.png', dpi=600,bbox_inches='tight')
            plt.close()
else:
    method = ''
    dorners = ['_TD1_','_TD2_','_TD3_', '_TD5_','_TD6_','_TD8_']
    # '_TD1_','_TD2_',
    # '_TD2_',
    for dorner in dorners: 
        work_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/bulk_dataset1/" + dorner +"/SVG/"
        print(dorner)
        
        if method == "spatalk":
            out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + dorner  + "/venn/spatalk/"
            dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/bulk_dataset1/' + dorner + '/Result_Repair/spatalk/'
        elif method == "sprout":
            out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + dorner + "/venn/sprout/"
            dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/bulk_dataset1/' + dorner + '/Result_Repair/sprout/'
        else:
            out_dir = "/home/xuezhengyang/data6/02-deconv_1/Script/Figure/" + dorner + "/venn/ST/"
            dir = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/filtered/st/st' + dorner + '.h5ad'
            
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            
        idata = sc.read(work_dir + 'meta_adata.h5ad')
        
        exp = idata.to_df()
        
        # if method == "spatalk" or method == "sprout":
        #     exp_path = dir + 'exp.csv'
        #     genes_path = dir + 'genes.csv'
            
        #     exp = pd.read_csv(exp_path)
        #     genes =pd.read_csv(genes_path)
            
        #     exp.index = genes.values[:,0]
        #     exp = exp.T
        # else:
        #     exp = sc.read(dir)
        #     exp.obs_names_make_unique()
        #     exp.var_names_make_unique()
        #     exp = exp.to_df()
            
        I_value = idata.uns['moranI']["I"].dropna()
        
        C_value = idata.uns['gearyC']["C"].dropna()
        
        filt = list(I_value.index)
        
        idata = idata[:,filt]
        
        # exp = exp.loc[idata.obs_names,idata.var_names]
        
        index = idata.var_names
        
        SOMDE_index = idata.uns['SOMDE'].index
        SOMDE_index = index[SOMDE_index]
        SOMDE_index = list(SOMDE_index)
        
        SPARKX_df = idata.uns['SPARKX']['adjustedPval']
        SPARKX_df = SPARKX_df[(SPARKX_df<threshold)]
        SPARKX_index = list(SPARKX_df.index)
        if len(SPARKX_index) == 0:
            SPARKX_df = idata.uns['SPARKX']['combinedPval']
            SPARKX_df = SPARKX_df[(SPARKX_df<threshold)]
            SPARKX_index = list(SPARKX_df.index)
            
        SpatialDE_df = idata.uns['SpatialDE']['padj']
        SpatialDE_df = SpatialDE_df[(SpatialDE_df<threshold)]
        SpatialDE_index = list(SpatialDE_df.index)
        
        nnSVG_index = idata.uns['nnSVG'].index
        # nnSVG_index = index[nnSVG_index]
        nnSVG_index = list(nnSVG_index)
        
        idata.uns['scGCO'].index = idata.uns['scGCO']["Unnamed: 0"]
        scGCO_df = idata.uns['scGCO']['p_value']
        values = idata.uns['scGCO']['p_value'].values
        for i in range(0,len(values)):
            values[i] = values[i][1:-1]
            values[i] = values[i].split(',')
            values[i] = float(values[i][0])
        idata.uns['scGCO']['p_value'] = values
        scGCO_df = scGCO_df[(scGCO_df<threshold)]
        scGCO_index = list(scGCO_df.index)
        # if scGCO_index.size == 0:
        #     scGCO_df = idata.uns['nnSVG']['combinedPval']
        #     scGCO_df = scGCO_df[(scGCO_df<threshold)]
        #     scGCO_index = list(scGCO_df.index)
        
        
        
        df = pd.DataFrame(index=index, columns=["SOMDE", "SPARKX", "SpatialDE","nnSVG","scGCO", "I", "C"])
        
        df[["SOMDE", "SPARKX", "SpatialDE","nnSVG","scGCO", "I", "C"]] = 0
        
        df.loc[SOMDE_index,"SOMDE"] = 1
        
        df.loc[SOMDE_index,"nnSVG"] = 1
        
        df.loc[SPARKX_index,"SPARKX"] = 1
        
        df.loc[SPARKX_index,"scGCO"] = 1
        
        df.loc[SpatialDE_index,"SpatialDE"] = 1
        
        df["C"] = C_value
        df["I"] = I_value
        
        corr = idata.obs[['row','col']]
        
        exp = pd.concat([corr,exp],axis=1)
        
        exp.index.name='id'
        df.index.name='id'
        
        exp.to_csv(out_dir + "exp.csv")
        df.to_csv(out_dir + "meta.csv")