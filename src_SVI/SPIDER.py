import pandas as pd
import numpy as np
import anndata 
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from mycolorpy import colorlist as mcp
from os.path import exists
from os import mkdir
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.preprocessing import MinMaxScaler
from pathlib import Path
PACKAGEDIR = Path(__file__).parent.absolute()


from . import svi
# from . import preprocess

class SPIDER():
    def __init__(self):
        self.svi = svi
        # self.pp = preprocess
        pass


    def prep(self,
            adata_input, work_dir, 
            no_spatalk=False,
            cluster_key='type', 
            pos_key=['row', 'col'], 
            is_human=True, 
            dist_offset=np.sqrt(2), 
            grid=True, 
            imputation=True,
            overwrite=False,
    ):
        adata = adata_input.copy()
        del adata_input
        # detect cci
        if not no_spatalk:
            out_f = f'{work_dir}/spatalk'
            if overwrite | (not exists(f'{out_f}_lrpair.csv')):
                preprocess.cci_spatalk(adata, work_dir, cluster_key, pos_key, is_human, out_f)
            else:
                print(f'{out_f}_lrpair.csv already exists, skipping spatalk.')
            lr_raw = pd.read_csv(f'{out_f}_lrpair.csv', index_col=0)
        else:
            lr_raw = preprocess.load_lr_df(is_human)     
        if imputation:
            'Running imputation with MAGIC'
            preprocess.impute_MAGIC(adata)

        # idata generation
        lr_df, adata = preprocess.subset(adata, lr_raw[['ligand', 'receptor']].drop_duplicates())
        pairs = preprocess.find_pairs(adata, pos_key, grid, dist_offset)
        pairs_meta = preprocess.meta(adata, pos_key, cluster_key, pairs)
        score = preprocess.score(adata, lr_df, pairs)
        idata = preprocess.idata_construct(score, pairs_meta, lr_df, lr_raw, adata)
        return idata

    def find_svi(self, idata, out_f, abstract=True, overwrite=False, n_neighbors=10, threshold=0.01):
        if not exists(out_f):
            print(f'Creating folder {out_f}')
            mkdir(out_f)
        if abstract:
            som, idata, meta_idata = svi.abstract(idata, n_neighbors)
            svi.find_svi(meta_idata,out_f,overwrite, som=som) #generating results
            svi.find_svi(idata, out_f,overwrite) #reading restuls
            svi_df, svi_df_strict = svi.combine_SVI(idata,threshold=threshold)
            svi_df.to_csv(f'{out_f}all.csv')
            if (overwrite) | (not exists(f'{out_f}pattern.csv')):
                svi.SVI_patterns(meta_idata, svi_df_strict)
                pd.DataFrame(meta_idata.obsm['pattern_score']).to_csv(f'{out_f}pattern.csv')
                meta_idata.var.to_csv(f'{out_f}membership.csv')
            else:
                svi.SVI_patterns(meta_idata, svi_df_strict)
                pd.DataFrame(meta_idata.obsm['pattern_score']).to_csv(f'{out_f}pattern.csv')
                # meta_idata.obsm['pattern_score'] = pd.read_csv(f'{out_f}pattern.csv', index_col=0).to_numpy()
                meta_idata.var.to_csv(f'{out_f}membership.csv')
                # meta_idata.var = pd.read_csv(f'{out_f}membership.csv', index_col=0)
            svi.meta_pattern_to_idata(idata, meta_idata)
            idata.var = meta_idata.var
            return idata, meta_idata
        else:
            svi.find_svi(idata, out_f,overwrite) #generating results
            svi_df, svi_df_strict = svi.combine_SVI(idata,threshold=threshold)
            svi_df.to_csv(f'{out_f}all.csv')
            if len(svi_df_strict) < 10:
                svi_df, svi_df_strict = svi.combine_SVI_Fisher(idata,threshold=threshold)
                print('Detected SVI number is less than 10, falling back to relaxed filtering.')
            if len(svi_df_strict) < 10:
                svi_df, svi_df_strict  = svi.combine_SVI_somde(idata,threshold=threshold)
                print('Detected SVI number is less than 10, falling back to use SOMDE result only.')
            if (overwrite) | (not exists(f'{out_f}pattern.csv')):
                svi.SVI_patterns(idata, svi_df_strict)
                pd.DataFrame(idata.obsm['pattern_score']).to_csv(f'{out_f}pattern.csv')
                idata.var.to_csv(f'{out_f}membership.csv')
            else:
                idata.obsm['pattern_score'] = pd.read_csv(f'{out_f}pattern.csv', index_col=0)
                idata.var = pd.read_csv(f'{out_f}membership.csv', index_col=0)                
            return idata
    
    def vis_pattern(self, idata, label='', pos_key=[], obsm_key='pattern_score', traj_coords=[], show_SVI=0):
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



            
    def interface_clustering_test(self, meta_idata, out_f):
        from pyseat.SEAT import SEAT
        import umap, umap.plot

        ns = [20, 50, 100]
        ds = [0.5, 0.8, 0.99]
        for n in ns:
            for d in ds:
                reduced = umap.UMAP(random_state=42, n_neighbors=n, min_dist=d, n_components=2).fit(meta_idata.obsm['patterns'])
                X = pd.DataFrame(reduced.embedding_, index=meta_idata.obs_names)
                x_std = np.std(X.to_numpy(), axis=None)
                print("in loop")
                for aff in ["gaussian_kernel", "linear_kernel", "cosine_similarity", "knn_neighbors_from_X", "laplacian_kernel"]:
                    for spas in ["knn_neighbors_from_X"]:
                        for n_nei in [20, 30, 50, 100]:
                            if (aff != 'linear_kernel') & ('kernel' in aff):
                                gm_scales = [0.01, 0.1, 1, 10, 100, 1000]
                            else:
                                gm_scales = [1]
                            
                            for gm_scale in gm_scales:
                                seat = SEAT(affinity=aff,
                                            sparsification=spas,
                                            objective="SE",
                                            kernel_gamma=gm_scale*x_std,
                                            n_neighbors=n_nei,
                                            max_k = 10,
                                            strategy="bottom_up")
                                seat_clusters = seat.fit_predict(X)
                                cluster_df = seat.ks_clusters
                                cluster_df.index = meta_idata.obs_names
                                cluster_df['Interface_Club'] = seat.clubs
                                cluster_df['Interface_Subpopulation'] = seat.labels_
                                cluster_df.to_csv(f'{out_f}/{"_".join("%s" %id for id in [n, d, aff, n_nei, gm_scale])}.csv')

    def interface_clustering(self, idata, meta_idata, n, d, aff, n_nei, gm_scale):
        from pyseat.SEAT import SEAT
        import umap, umap.plot
        
        reduced = umap.UMAP(random_state=42, n_neighbors=n, min_dist=d, n_components=2).fit(meta_idata.obsm['patterns'])
        X = pd.DataFrame(reduced.embedding_, index=meta_idata.obs_names)
        x_std = np.std(X.to_numpy(), axis=None)

        seat = SEAT(affinity=aff,
                    sparsification="knn_neighbors_from_X",
                    objective="SE",
                    n_neighbors=n_nei,
                    kernel_gamma=gm_scale*x_std,
                    max_k=10,
                    strategy="bottom_up")
        seat.fit_predict(X)
        
        cluster_df = seat.ks_clusters
        cluster_df.index = meta_idata.obs_names
        cluster_df['Interface_Club'] = seat.clubs
        cluster_df['Interface_Subpopulation'] = seat.labels_
        
        idata.obs[cluster_df.columns] = cluster_df.iloc[idata.obs['som_node'].to_numpy()].to_numpy()
        return reduced, seat, idata
    
    def labeled_spot_interface(self, idata, label_key='Interface_Subpopulation'):
        belonging = {}
        cells = idata.uns['cell_meta'].index
        for i in cells:
            belonging[i] = []
        for pair in idata.obs.reset_index()[[label_key,'A', 'B']].to_numpy():
            belonging[pair[1]].append(pair[0])
            belonging[pair[2]].append(pair[0])
        import collections
        df = pd.concat([pd.Series(collections.Counter(belonging[c])) for c in cells], axis=1).T
        df = df.fillna(0)
        df.index = cells
        idata.uns['interface_set'] = df
        return df
    
    def scored_spot_interface(self, idata):       
        belonging = {}
        cells = idata.uns['cell_meta'].index
        for i in cells:
            belonging[i] = []
        for pair in idata.obs.reset_index()[['index','A', 'B']].to_numpy():
            belonging[pair[1]].append(pair[0])
            belonging[pair[2]].append(pair[0])
        score = pd.DataFrame(idata.obsm['pattern_score'], index=idata.obs_names)
        df = pd.concat([score.loc[belonging[c]].mean() for c in cells], axis=1).T     
        df.index = cells
        idata.uns['cell_pattern'] = df
        print(f'Added key cell_pattern in idata.uns')
        return df

    def interaction_spot_interface(self, idata):       
        belonging = {}
        cells = idata.uns['cell_meta'].index
        for i in cells:
            belonging[i] = []
        for pair in idata.obs.reset_index()[['index','A', 'B']].to_numpy():
            belonging[pair[1]].append(pair[0])
            belonging[pair[2]].append(pair[0])
        score = idata.to_df()
        df = pd.concat([score.loc[belonging[c]].mean() for c in cells], axis=1).T     
        df.index = cells
        idata.uns['cell_pattern'] = df
        print(f'Added key cell_interaction in idata.uns')
        return df
    
    
    def viz_interface_pattern(self, idata, pos_key=['row', 'col'], label='',  histology=None):
        plt.figure(figsize=(15, 15))
        if label != '':
            plt.subplot(6, 5, 1)
            plt.scatter(idata.uns['cell_meta'][pos_key[0]], idata.uns['cell_meta'][pos_key[1]], c=idata.uns['cell_meta'][label].astype("category").cat.codes, s=1)
            plt.axis('equal')
            base = 2
        else:
            base = 1
        for i in range(idata.uns['cell_pattern'].shape[1]):
            plt.subplot(6, 5, i + base)
            plt.scatter(idata.uns['cell_meta'][pos_key[0]], idata.uns['cell_meta'][pos_key[1]], c=idata.uns['cell_pattern'][i], s=1)
            plt.axis('equal')
            if histology:
                plt.title('Pattern {} - {} LRIs'.format(i, histology.query('pattern == @i').shape[0] ))
            plt.colorbar(ticks=[])

    def viz_interface_pattern_on_spot(self, idata, pattern_id, pos_key=['row', 'col'], label='', ax=None):
        arr = idata.obs[['A_row', 'A_col', 'B_row', 'B_col']]
        arr['color'] = mcp.gen_color_normalized(cmap="RdBu_r",data_arr=idata.obsm['pattern_score'][:, pattern_id])
        arr = arr.to_numpy()
        for j in arr:
            ax.plot([j[0], j[2]], [j[1], j[3]], color = j[4], alpha=0.8, linewidth=3)
        # draw cells
        fig = ax.scatter(x=idata.uns['cell_meta'][pos_key[0]], y=idata.uns['cell_meta'][pos_key[1]], c=idata.uns['cell_meta'][label].astype('category').cat.codes, zorder=2, s=1, cmap='tab10')
        ax.axis('equal')
        ax.axis('off')
        ax.set_title(f'Pattern {pattern_id}')
        return fig

    def viz_LRI_on_spot(self, idata, id_LRI, pos_key=['row', 'col'], label='', ax=None):
        arr = idata.obs[['A_row', 'A_col', 'B_row', 'B_col']]
        arr['color'] = mcp.gen_color_normalized(cmap="RdBu_r",data_arr=idata.to_df()[id_LRI])
        arr = arr.to_numpy()
        for j in arr:
            ax.plot([j[0], j[2]], [j[1], j[3]], color = j[4], alpha=0.8, linewidth=3)
        # draw cells
        fig = ax.scatter(x=idata.uns['cell_meta'][pos_key[0]], y=idata.uns['cell_meta'][pos_key[1]], c=idata.uns['cell_meta'][label].astype('category').cat.codes, zorder=2, s=1, cmap='tab10')
        ax.axis('equal')
        ax.axis('off')
        ax.set_title(id_LRI)
        return fig

    
    def spot_clust(self, idata, adata, label, portion=0.1, ideal_k=None, random_state=0):
        from pyseat.SEAT import SEAT
        import umap, umap.plot
        
        masked_target = idata.uns['cell_meta'][label].copy()
        np.random.seed(random_state)
        masked_target[np.random.choice(len(idata.uns['cell_meta'][label]), size=int(len(idata.uns['cell_meta'][label]) * portion), replace=False)] = -1
        embedding = umap.UMAP(n_neighbors=100, min_dist=1).fit_transform(idata.uns['cell_pattern'], y=masked_target)
        adata.obsm['X_umap'] = embedding
        sc.pp.neighbors(adata, n_neighbors=20, use_rep='X_umap')
        sc.tl.draw_graph(adata)
        flag = 1
        res_max = 10
        res_min = 0.0001
        res = 0.05
        if not ideal_k:
            k = len(idata.uns['cell_meta']['layer_guess'].unique())
        while (flag):
            sc.tl.leiden(adata, resolution=res)
            print(len(adata.obs['leiden'].unique()), res)
            if (len(adata.obs['leiden'].unique()) == k) | (res < 0) | (res > 10):
                flag = 0
            elif len(adata.obs['leiden'].unique()) < k:
                res_min = res
                res = (res+res_max)/2
            else:
                res_max = res
                res = (res+res_min)/2
        return adata.obs['leiden'], embedding


    def relabel_interface(self, idata, clust_key):
        node_labels_text = idata.uns['cell_meta'][clust_key]
        idata.obs[f'A_label_{clust_key}'] = node_labels_text.loc[idata.obs['A']].astype(str).to_numpy()
        idata.obs[f'B_label_{clust_key}'] = node_labels_text.loc[idata.obs['B']].astype(str).to_numpy()
        node_labels = idata.uns['cell_meta'][clust_key].astype('category').cat.codes
        idata.obs[f'A_label_int_{clust_key}'] = node_labels.loc[idata.obs['A']].to_numpy()
        idata.obs[f'B_label_int_{clust_key}'] = node_labels.loc[idata.obs['B']].to_numpy()
        idata.obs['label_1'] = idata.obs[f'A_label_int_{clust_key}'].astype(str) + idata.obs[f'B_label_int_{clust_key}'].astype(str)
        idata.obs['label_2'] = idata.obs[f'B_label_int_{clust_key}'].astype(str) + idata.obs[f'A_label_int_{clust_key}'].astype(str)
        idata.obs['label_int_{clust_key}'] = idata.obs[['label_1', 'label_2']].astype(int).max(axis=1).astype(str).astype('category')
        label_1 = idata.obs[f'A_label_{clust_key}'].astype(str) + '_' + idata.obs[f'B_label_{clust_key}'].astype(str).to_numpy()
        label_2 = idata.obs[f'B_label_{clust_key}'].astype(str) + '_' + idata.obs[f'A_label_{clust_key}'].astype(str).to_numpy()
        pick = idata.obs[['label_1', 'label_2']].astype(int).idxmax(axis=1).to_numpy()
        text_label = [label_1[i] if x=='label_1' else label_2[i] for i,x in enumerate(pick)]
        idata.obs[f'label_{clust_key}'] = text_label
        idata.obs[f'label_{clust_key}'] = idata.obs[f'label_{clust_key}'].astype('category')
        print(f'Added key label_{clust_key} in idata.obs')
    
    def list_pattern_LRI(self, histology_results, short=True):
        for i in histology_results.sort_values('pattern').pattern.unique():
            print('Pattern {}'.format(i))
            if short:
                print(', '.join(histology_results.query('pattern == @i').sort_values('membership')['g'][:10].tolist()))
            else:
                print(', '.join(histology_results.query('pattern == @i').sort_values('membership')['g'].tolist()))
            print()
            
    def labeled_spot_interface(self, idata, label_key='Interface_Subpopulation', normalization=True):
        belonging = {}
        cells = idata.uns['cell_meta'].index
        for i in cells:
            belonging[i] = []

        for pair in idata.obs.reset_index()[[label_key,'A', 'B']].to_numpy():
            belonging[pair[1]].append(pair[0])
            belonging[pair[2]].append(pair[0])
            
        import collections
        df = pd.concat([pd.Series(collections.Counter(belonging[c])) for c in cells], axis=1).T
        df = df.fillna(0)
        df.index = cells
        idata.uns['interface_set'] = df
        if normalization:
            df_nor = df.div(df.sum(axis=1), axis=0)
            idata.uns['interface_set_normalized'] = df_nor
            return df_nor
        else:
            return df

    def load_pathway(self, is_human=True):
        pw = pd.read_csv(f'{PACKAGEDIR}/lr_pair_data/pathways.tsv', sep='\t', index_col=0)
        species = 'Human' if is_human else 'Mouse'
        pw = pw[pw.species==species].drop_duplicates()
        pw.index = pw.src + '_' + pw.dest
        return pw
        

    def pathway_annotation(self, idata, is_human=True):
        pw = self.load_pathway(is_human)
        pw = pw.loc[pw.index.isin(idata.var_names)]
        df = pd.get_dummies(pw.pathway)
        df = df.groupby(df.index).sum()
        idata.varm['pathway'] = pd.DataFrame(0, columns=df.columns, index=idata.var_names)
        idata.varm['pathway'].loc[df.index] = df
        print(f'Added key pathway in idata.varm')

        
    def pathway_annotation_list(self, LRI_list, is_human=True):
        pw = self.load_pathway(is_human)
        pw = pw.loc[pw.index.isin(LRI_list)]
        df = pd.get_dummies(pw.pathway)
        df = df.groupby(df.index).sum()
        return df

    def pathway_prep(self, idata, is_human=True):
        lr_raw = preprocess.load_lr_df(is_human=is_human)
        lr_raw.index = lr_raw.ligand + '_' + lr_raw.receptor
        pw = self.load_pathway(is_human)
        pw.index = pw.src + '_' + pw.dest
        pw = pw.loc[np.intersect1d(lr_raw.index, pw.index)]
        custom = {}
        for ptw in pw.pathway.unique():
            custom[ptw] = pw[pw.pathway == ptw].index.tolist()
        background = np.intersect1d(idata.var_names, pw.index).tolist()
        return custom, background
    
    def pathway_prep_custom_background(self, background, is_human=True):
        lr_raw = preprocess.load_lr_df(is_human=is_human)
        lr_raw.index = lr_raw.ligand + '_' + lr_raw.receptor
        pw = self.load_pathway(is_human)
        pw.index = pw.src + '_' + pw.dest
        pw = pw.loc[np.intersect1d(lr_raw.index, pw.index)]
        custom = {}
        for ptw in pw.pathway.unique():
            custom[ptw] = pw[pw.pathway == ptw].index.tolist()
        background = np.intersect1d(background, pw.index).tolist()
        return custom, background

    def enrichment_v2(self, idata, is_human=True, groupby='label', subset=[], order=[], custom=None, background=None):
        import gseapy
        if not (custom and background):
            custom, background = self.pathway_prep(idata, is_human=is_human)
        var = idata.var
        if len(subset) != 0:
            var = var[var[groupby].isin(subset)]
        arr = []
        for i in var.sort_values(groupby)[groupby].unique():
            lri_list = var[var[groupby] ==i].index.to_numpy().tolist()
            arr.append([str(i), lri_list])
        dfs = []
        for i in range(len(arr)):
            try:
                enr_res = gseapy.enrichr(gene_list=arr[i][1],
                                gene_sets=custom,
                                background=background)
                                # cutoff = cutoff)
                enr_res.res2d[groupby] = arr[i][0]
                dfs.append(enr_res.res2d)
            except Exception as e:
                print(e)
                continue
        merged_df = pd.concat(dfs)
        if len(order)!=0:
            from string import ascii_lowercase
            rename_dict = {}
            count = 0
            for x in order:
                rename_dict[str(x)] = f'({ascii_lowercase[count]}) {x}'
                count += 1
            merged_df[f'ordered_{groupby}'] = merged_df[groupby]
            merged_df[f'ordered_{groupby}'] = merged_df[groupby].map(rename_dict)
        else:
            merged_df[f'ordered_{groupby}'] = merged_df[groupby]
        return merged_df, arr

    def enrichment_gene(self, idata, is_human=True, groupby='label', subset=[], order=[], custom=None, background=None, custom_pathwaydb=[]):
        import gseapy
        var = idata.var
        if len(subset) != 0:
            var = var[var[groupby].isin(subset)]
        arr = []
        for i in var.sort_values(groupby)[groupby].unique():
            lri_list = var[var[groupby] ==i].index.to_numpy().tolist()
            gene_list = np.unique(np.concatenate([x.split('_') for x in lri_list])).tolist()
            arr.append([str(i), gene_list])
        if len(custom_pathwaydb) == 0:
            pathway_db = ['KEGG_2021_Human' if is_human else 'KEGG_2019_Mouse']
        else:
            pathway_db = custom_pathwaydb
        organism = 'Human' if is_human else 'Mouse'
        dfs = []
        for i in range(len(arr)):
            try:
                enr_res = gseapy.enrichr(gene_list=arr[i][1],
                                organism=organism,
                                gene_sets=pathway_db)
                enr_res.res2d[groupby] = arr[i][0]
                dfs.append(enr_res.res2d)
            except Exception as e:
                print(e)
                continue
        merged_df = pd.concat(dfs)
        if len(order)!=0:
            from string import ascii_lowercase
            rename_dict = {}
            count = 0
            for x in order:
                rename_dict[str(x)] = f'({ascii_lowercase[count]}) {x}'
                count += 1
            merged_df[f'ordered_{groupby}'] = merged_df[groupby]
            merged_df[f'ordered_{groupby}'] = merged_df[groupby].map(rename_dict)
        else:
            merged_df[f'ordered_{groupby}'] = merged_df[groupby]
        return merged_df, arr
        
        

    def enrichment(self, custom, background, histology_results, groupby='pattern', cutoff=1, top=None, order=[], group=[]):
        import gseapy
        dfs=[]
        histology_results = histology_results[histology_results.g.isin(np.intersect1d(histology_results.g, background))]
        arr = []
        if len(group) == 0:
            for i in histology_results.sort_values(groupby).pattern.unique():
                lri_list = histology_results[histology_results[groupby] ==i].sort_values('membership', ascending=False)['g'].to_numpy()
                if top:
                    lri_list = lri_list[:top]
                lri_list = lri_list.tolist()
                arr.append([str(i), lri_list])
        else:
            count = 0
            for x in order:
                index_m  = np.where(group == x)[0].tolist()
                lri_lists = []
                for i in index_m:
                    lri_lists.append(histology_results.query('@groupby == @i').sort_values('membership', ascending=False)['g'].tolist())
                lri_list = np.concatenate(lri_lists)
                if top:
                    lri_list = lri_list[:top]
                lri_list = lri_list.tolist()
                arr.append([f'Merged groupby {count}: {index_m}', lri_list])
                count += 1

        for i in range(len(arr)):
            try:
                enr_res = gseapy.enrichr(gene_list=arr[i][1],
                                gene_sets=custom,
                                background=background,
                                cutoff = cutoff)
                enr_res.res2d['group'] = arr[i][0]
                dfs.append(enr_res.res2d)
            except:
                print(i)
                continue
        merged_df = pd.concat(dfs)
        if len(order)!=0 and len(group) == 0:
            from string import ascii_lowercase
            rename_dict = {}
            count = 0
            for x in order:
                rename_dict[str(x)] = f'({ascii_lowercase[count]}) {x}'
                count += 1
            merged_df['ordered_group'] = merged_df.group
            merged_df['ordered_group'] = merged_df.group.map(rename_dict)
        else:
            merged_df['ordered_group'] = merged_df.group
        return merged_df, arr

    def vis_enrichment(self, df, x_key='ordered_group', figsize=(5, 5), size=3, save=None, cutoff=0.1, top_term=10):
        import gseapy
        ax = gseapy.dotplot(df,figsize=figsize, 
                        x=x_key, 
                        top_term=top_term,
                        column='P-value',
                        cutoff=cutoff,
                        cmap = plt.cm.plasma, 
                        xticklabels_rot=0,
                        size=size, 
                        ofname=save,
                        show_ring=True)
        if not save:
            # ax.set_xlabel("")
            plt.show()

    def add_diff(self, idata, list1, list2, key_to_add):
        scaler = MinMaxScaler()
        exp1=scaler.fit_transform(np.sqrt(idata[:, idata.var_names.isin(list1)].to_df()).sum(axis=1).to_numpy().reshape(-1, 1)).flatten()
        exp2=scaler.fit_transform(np.sqrt(idata[:, idata.var_names.isin(list2)].to_df()).sum(axis=1).to_numpy().reshape(-1, 1)).flatten()
        idata.obs[key_to_add]=exp1-exp2
        print(f'Added key {key_to_add} in idata.uns')

        


    # def spot_clust_old(self, idata, label, interface_weight=0.5, n=20, d=0.5, aff='gaussian_kernel', n_nei=20, gm_scale=1, ideal_k=None):
    #     from pyseat.SEAT import SEAT
    #     import umap, umap.plot
    #     method_label_df = pd.get_dummies(idata.uns['cell_meta'][label])
    #     X = pd.concat([method_label_df, idata.uns['cell_pattern']], axis=1)
    #     w = np.concatenate((np.repeat( (1-interface_weight)/method_label_df.shape[1], method_label_df.shape[1]), 
    #                     np.repeat( interface_weight/idata.uns['cell_pattern'].shape[1], idata.uns['cell_pattern'].shape[1])))
    #     A = X.to_numpy()@np.diag(w)
    #     reduced = umap.UMAP(random_state=42, n_neighbors=n, min_dist=d, n_components=2).fit(A)
    #     if not ideal_k:
    #         ideal_k = method_label_df.shape[1]
    #     seat = SEAT(affinity=aff,
    #         sparsification="knn_neighbors_from_X",
    #         objective="SE",
    #         n_neighbors=n_nei,
    #         max_k = ideal_k,
    #         strategy="bottom_up")
    #     seat.fit_predict(reduced.embedding_)
    #     return seat.ks_clusters.to_numpy()[:, -1], reduced
    
    # def spot_clust(self, idata, label, ideal_k=None, random_state=0):
    #     from pyseat.SEAT import SEAT
    #     import umap, umap.plot
        
    #     X_scaled = MinMaxScaler(feature_range=(-1, 1)).fit_transform(idata.uns['cell_pattern'])
    #     masked_target = idata.uns['cell_meta'][label].copy().astype(np.int8)
    #     np.random.seed(random_state) 
    #     masked_target[np.random.choice(len(idata.uns['cell_meta'][label]), size=int(len(idata.uns['cell_meta'][label]) * 0.2), replace=False)] = -1
    #     reduced = umap.UMAP(random_state=42, n_neighbors=50, min_dist=0.1, n_components=2).fit(X_scaled, y=masked_target)
    #     if not ideal_k:
    #         ideal_k = len(idata.uns['cell_meta'][label].unique())
    #     seat = SEAT(affinity="gaussian_kernel",
    #         sparsification="knn_neighbors_from_X",
    #         objective="SE",
    #         n_neighbors=15,
    #         max_k = ideal_k,
    #         strategy="bottom_up")
    #     seat.fit_predict(reduced.embedding_)
    #     return seat.ks_clusters.to_numpy()[:, -1], reduced

