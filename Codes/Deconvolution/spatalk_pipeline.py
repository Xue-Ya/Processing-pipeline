import pandas as pd
import os

def cci_spatalk(adata, work_dir, cluster_key, pos_key, is_human, out_f):
    count_f = f'{work_dir}/adata_count.csv'
    meta_f = f'{work_dir}/adata_meta.csv'
    
    df = adata.to_df()
    df.index = "C"+df.index
    df.to_csv(count_f)
    meta = adata.obs[pos_key+[cluster_key]].reset_index()
    meta.columns = ['cell', 'x', 'y', 'celltype']
    meta.cell = "C"+meta.cell
    if not pd.api.types.is_string_dtype(meta.celltype.dtype):
        meta.celltype = "T"+meta.celltype.astype('str')
    meta.to_csv(meta_f)
    species = 'Human' if is_human else 'Mouse'
    os.system(str(f'/bin/bash -c "source /etc/profile;module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Anaconda3/2022.05 R-bundle-Bioconductor/3.15-R-4.2.0;R -f ../src/run_spatalk.R {count_f} {meta_f} {species} {out_f}"'))