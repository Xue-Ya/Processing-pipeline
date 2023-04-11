#!/usr/bin/env /home/xuezhengyang/.conda/envs/FYP/bin/python
import scanpy as sc
import pandas as pd
import seaborn as sns

import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

sc.settings.verbosity = 1

# total_before = 'total6'
# total_now = 'total7'
# adata_path = '/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset3/niche_raw_filtered.h5ad'

datas = [
    "/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset2/sc/dataset_2_filtered.h5ad",
    "/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/sc_filtered.h5ad",
    "/home/xuezhengyang/data6/02-deconv_1/Script/Data/bulk_dataset1/sc_filtered.h5ad",
    "/home/xuezhengyang/data6/02-deconv_1/Script/notebook/check_Lin/dataset_3/dataset_3_filtered.h5ad",
    "/home/xuezhengyang/data6/02-deconv_1/Script/notebook/check_Lin/dataset_6/dataset_6_filtered.h5ad",
    "/home/xuezhengyang/data6/02-deconv_1/Script/notebook/check_Lin/dataset_8/dataset_8_filtered.h5ad",
    "/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset3/niche_raw_filtered.h5ad"
]

# adata_ref = sc.read( '/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/' + total_before + '.h5ad')
adata_ref = sc.read( '/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/total_rna_droped.h5ad')
total_now = 0

for adata_path in datas:
    
    logger.info(adata_ref.shape)
    # logger.info()
    # index = []
    # for i in range(0,len(list(adata_ref.obs_names))) :
    #     index.append(adata_ref.obs_names[i][0:-4])
    # adata_ref.obs_names = index
    # logger.info(adata_ref.obs_names)
    adata = sc.read(adata_path)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    var_names = adata_ref.var_names.intersection(adata.var_names)
    adata_ref = adata_ref[:, var_names]
    adata = adata[:, var_names]

    logger.info(var_names.shape)

    sc.pl.umap(adata_ref, color='ann_finest_level',save = '/insert/' + str(total_now) + '_before.png',show=False)

    sc.pp.pca(adata)
    sc.pp.neighbors(adata)

    sc.pp.pca(adata_ref)
    sc.pp.neighbors(adata_ref)

    adata_ref.uns['umap'] = {'params': {'a': 0.583030019901822, 'b': 1.3341669931033755}}

    sc.tl.ingest(adata, adata_ref, obs='ann_finest_level',inplace=True)

    sc.pl.umap(adata, color='ann_finest_level', save = '/insert/' + str(total_now) + '_adata.png',show=False)

    adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])

    adata_concat.obs.ann_finest_level = adata_concat.obs.ann_finest_level.astype('category')
    adata_concat.obs.ann_finest_level.cat.reorder_categories(adata_ref.obs.ann_finest_level.cat.categories, inplace=True)  # fix category ordering
    # adata_concat.uns['louvain_colors'] = adata_ref.uns['louvain_colors']  # fix category colors



    sc.pl.umap(adata_concat, color='ann_finest_level',save = '/insert/' + str(total_now) + '_after.png',show=False)

    index = list(adata_concat.obs_names)

    logger.info(adata_concat.obs_names[:200])
    
    for i in range(0,len(index)) :
        index[i] = index[i][0:-4]
    adata_concat.obs_names = index
    
    logger.info(adata_concat.shape)
    logger.info(adata_concat.obs_names[:200])
    
    adata_ref = adata_concat
    total_now = total_now + 1
    
    adata_concat.write('/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/' + 'res' + '.h5ad')

