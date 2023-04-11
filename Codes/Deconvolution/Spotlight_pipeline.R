
library(Seurat)
library(SPOTlight)

args = commandArgs()

if (length(args)==0) {
  stop("not enough input", call.=FALSE)
}

st_data <-  args[4]
st_meta <-  args[5]
sc_data <- args[6]
sc_celltype <- args[7]
# species <- args[8]
# recons_out_f <- args[9]
decon_out_f <- args[8]
# deconv_methoid <- args[11]

ref_data<- Seurat::CreateSeuratObject(sc_data)
ref_data<- Seurat::SCTransform(ref_data, verbose = F)
Seurat::Idents(ref_data)<- sc_celltype$celltype
cluster_markers_all <- Seurat::FindAllMarkers(object = ref_data, assay = "SCT", slot = "data", verbose = F, only.pos = TRUE)
ref_data$celltype <- sc_celltype$celltype
spotlight_ls <- SPOTlight::spotlight_deconvolution(se_sc = ref_data, counts_spatial = as.matrix(st_data), clust_vr = "celltype", cluster_markers = cluster_markers_all)
spotlight_ls<- as.data.frame(spotlight_ls[[2]])
spotlight_ls<- spotlight_ls[,-ncol(spotlight_ls)]
st_coef <- as.matrix(st_coef)

write.csv(st_coef ,paste0(decon_out_f,'/RCTD.csv',sep = "",collapse = NULL),row.names = FALSE)