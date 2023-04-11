import sys
import os

class Deconvolutions:
    def __init__(self, sc_dir = None, RNA_h5ad = None, RNA_h5Seurat = None, Spatial_file = None, Spatial_h5ad = None, Spatial_h5Seurat = None, celltype_key = None, celltype_file = None, my_python_path = None, output_path = None,work_dir= None):
        """
            @author: wen zhang
            This function integrates spatial and scRNA-seq data to predictes the celltype deconvolution of the spots.
            
            A minimal example usage:
            Assume we have (1) scRNA-seq data file named RNA_h5ad or RNA_h5Seurat
            (2) spatial transcriptomics data file named Spatial_h5ad or Spatial_h5Seurat
            (3) celltype annotataion title in scRNA-seq data file
            
            >>> import Benchmarking.DeconvolutionSpot as DeconvolutionSpot
            >>> test = DeconvolutionSpot.Deconvolutions(RNA_file, RNA_h5ad, RNA_h5Seurat, Spatial_file, Spatial_h5ad, Spatial_h5Seurat, celltype_key, celltype_file, my_python_path, output_path)
            >>> Methods = ['Cell2location','SpatialDWLS','RCTD','STRIDE','Stereoscope','Tangram','DestVI', 'Seurat', 'SPOTlight', 'DSTG']
            >>> Result = test.Dencon(Methods)
            
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
        
        self.sc_dir = sc_dir
        self.RNA_h5ad = RNA_h5ad
        self.RNA_h5Seurat = RNA_h5Seurat
        self.Spatial_file = Spatial_file
        self.Spatial_h5ad = Spatial_h5ad
        self.Spatial_h5Seurat = Spatial_h5Seurat
        self.celltype_key = celltype_key
        self.celltype_file = celltype_file
        self.my_python_path = my_python_path
        self.output_path = output_path
        self.work_dir= work_dir
    

    
    def Dencon(self, need_tools):
        if "Tangram" in need_tools:
            RNA_h5ad = self.RNA_h5ad
            Spatial_h5ad = self.Spatial_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python Codes/Deconvolution/Tangram_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)

        if "DestVI" in need_tools:
            RNA_h5ad = self.RNA_h5ad
            Spatial_h5ad = self.Spatial_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python Codes/Deconvolution/DestVI_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)
            
        # if "RCTD" in need_tools:
        #     count_f = self.work_dir + '/Spatalk/spatalk_st_gene.csv'
        #     meta_f = self.work_dir + '/Spatalk/spatalk_st_corr.csv'
        #     sc_count = self.sc_dir + '/Spatalk/spatalk_sc_gene_sub.csv'
        #     sc_celltypes = self.sc_dir + '/Spatalk/spatalk_sc_celltype_sub.csv'
        #     output_path = self.output_path
            
        #     os.system(str(f'/bin/bash -c "source /etc/profile;module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Anaconda3/2022.05 R-bundle-Bioconductor/3.15-R-4.2.0;R -f /home/xuezhengyang/data6/02-deconv_1/Script/Codes/Deconvolution/RCTD_pipeline.R {count_f} {meta_f} {sc_count} {sc_celltypes} {output_path}"'))
            
        # if "SPOTlight" in need_tools:
        #     count_f = self.work_dir + '/Spatalk/spatalk_st_gene.csv'
        #     meta_f = self.work_dir + '/Spatalk/spatalk_st_corr.csv'
        #     sc_count = self.sc_dir + '/Spatalk/spatalk_sc_gene_sub.csv'
        #     sc_celltypes = self.sc_dir + '/Spatalk/spatalk_sc_celltype_sub.csv'
        #     output_path = self.output_path
            
        #     os.system(str(f'/bin/bash -c "source /etc/profile;module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Anaconda3/2022.05 R-bundle-Bioconductor/3.15-R-4.2.0;R -f /home/xuezhengyang/data6/02-deconv_1/Script/Codes/Deconvolution/Spotlight_pipeline.R {count_f} {meta_f} {sc_count} {sc_celltypes} {output_path}"'))
            
        # if "deconvSeq" in need_tools:
        #     RNA_h5ad = self.RNA_h5ad
        #     Spatial_h5ad = self.Spatial_h5ad
        #     celltype_key = self.celltype_key
        #     output_path = self.output_path
        #     os.system('python Codes/Deconvolution/deconvSeq_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)
        
        # if "stereoscope" in need_tools:
        #     RNA_h5ad = self.RNA_h5ad
        #     Spatial_h5ad = self.Spatial_h5ad
        #     celltype_key = self.celltype_key
        #     output_path = self.output_path
        #     os.system('python Codes/Deconvolution/stereoscope_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)
        
        # if "cell2location" in need_tools:
        #     RNA_h5ad = self.RNA_h5ad
        #     Spatial_h5ad = self.Spatial_h5ad
        #     celltype_key = self.celltype_key
        #     output_path = self.output_path
        #     os.system('python Codes/Deconvolution/cell2location_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)
