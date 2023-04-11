import sys
import os
import pandas as pd
import scanpy as sc

class Repair:
    def __init__(self, sc_dir = None, work_dir = None, spaTalk_dir = None,conv = None, conv_path = None, celltype_key = None,  output_path = None):
        """           
            Parameters
            -------
            
            RNA_file : str
            scRNA-seq data count file.
            
            RNA_h5ad : str
            scRNA-seq data file with h5ad format.
            
            RNA_h5Seurat : str
            scRNA-seq data file with h5Seurat format.
            
            Spatial_file : str
            Spatial data count file.
            
            Spatial_h5ad : str
            Spatial data file with h5ad format.
            
            Spatial_h5Seurat : str
            Spatial data file with h5Seurat format.
            
            celltype_key : str
            celltype annotataion title in scRNA-seq data h5ad file or h5Seurat file
            
            celltype_file : str
            celltype annotataion file
            
            my_python_path : str
            which python path used for Cell2location
            
            output_path : str
            Outfile path
            
            """
        
        self.work_dir = work_dir
        self.sc_dir = sc_dir
        self.spaTalk_dir = spaTalk_dir
        self.celltype_key = celltype_key
        self.output_path = output_path
        self.conv_path = conv_path
        self.conv = conv
       
        
    def Repair(self, need_tools):
        if "sprout" in need_tools:
            # RNA_file = self.work_dir + '/sc_gene.csv'
            RNA_file = self.sc_dir + '/sc_gene_sub.csv'
            celltype_file = self.sc_dir + '/sc_meta_sub.csv'
        
            # RNA_file = self.sc_dir + '/sc_gene_healthy.csv'
            # celltype_file = self.sc_dir + '/sc_meta_healthy.csv'
            
            # RNA_file = self.sc_dir + '/sc_gene_cancer.csv'
            # celltype_file = self.sc_dir + '/sc_celltype_cancer.csv'
            
            Spatial_file = self.work_dir + '/st_gene.csv'
            # celltype_file = self.work_dir + '/sc_meta.csv'

            Spitial_meta = self.work_dir + '/st_corr.csv'
            output_path = self.output_path
            lr_file_path = '/home/xuezhengyang/data6/02-deconv_1/Script/SPROUT-main/LR/human_LR_pairs.txt'
            conv_path = self.conv_path
            
            os.system('python SPROUT-SPROUT_fast/sprout.py '  + Spatial_file+ ' '  + Spitial_meta+ ' ' +conv_path+ ' '  + RNA_file+ ' '  + celltype_file+ ' '  +lr_file_path + ' '+ output_path)
            
        if "spatalk" in need_tools:
            is_human = True
            out_f = self.output_path
            species = 'Human' if is_human else 'Mouse'
            count_f = self.work_dir + '/Spatalk/spatalk_st_gene.csv'
            meta_f = self.work_dir + '/Spatalk/spatalk_st_corr.csv'
            sc_count = self.sc_dir + '/Spatalk/spatalk_sc_gene_sub.csv'
            sc_celltypes = self.sc_dir + '/Spatalk/spatalk_sc_celltype_sub.csv'
            conv_path = self.conv_path
            conv = self.conv
            
            print(conv)
            
            if conv == 0:
            #This is for finding LR pairs
            # os.system(str(f'/bin/bash -c "source /etc/profile;module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Anaconda3/2022.05 R-bundle-Bioconductor/3.15-R-4.2.0;R -f /home/xuezhengyang/data6/02-deconv_1/Script/Codes/Reconstruct/run_spatalk_LR.R {count_f} {meta_f} {sc_count} {sc_celltypes} {species} {out_f}"'))
                os.system(str(f'/bin/bash -c "source /etc/profile;module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Anaconda3/2022.05 R-bundle-Bioconductor/3.15-R-4.2.0;R -f /home/xuezhengyang/data6/02-deconv_1/Script/Codes/Reconstruct/run_spatalk.R {count_f} {meta_f} {sc_count} {sc_celltypes} {species} {out_f} {conv_path}"'))
            else:
                sc_count = self.sc_dir + '/Spatalk/spatalk_sc_gene_healthy.csv'
                sc_celltypes = self.sc_dir + '/Spatalk/spatalk_sc_celltype_healthy.csv'
                
                # sc_count = self.sc_dir + '/Spatalk/spatalk_sc_gene_cancer.csv'
                # sc_celltypes = self.sc_dir + '/Spatalk/spatalk_sc_celltype_cancer.csv'
                
                os.system(str(f'/bin/bash -c "source /etc/profile;module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Anaconda3/2022.05 R-bundle-Bioconductor/3.15-R-4.2.0;R -f /home/xuezhengyang/data6/02-deconv_1/Script/Codes/Reconstruct/run_spatalk_deconv.R {count_f} {meta_f} {sc_count} {sc_celltypes} {species} {out_f} {conv_path} {conv}"'))
                