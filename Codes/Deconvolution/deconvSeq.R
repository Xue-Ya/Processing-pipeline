library(deconvSeq)

# st_data <-  args[4]
# st_meta <-  args[5]
# sc_data <- args[6]
# sc_celltype <- args[7]
# # species <- args[8]
# # recons_out_f <- args[9]
# decon_out_f <- args[8]
# # deconv_methoid <- args[11]

st_data <-  "/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/st_12_2/Spatalk/spatalk_st_gene.csv"
st_meta <-  "/home/xuezhengyang/data6/02-deconv_1/Script/Data/Dataset1/apart/st_12_2/Spatalk/spatalk_st_corr.csv"
sc_data <- "/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/Spatalk/spatalk_sc_gene_sub.csv"
sc_celltype <- "/home/xuezhengyang/data6/02-deconv_1/Script/Data/SC/Spatalk/spatalk_sc_celltype_sub.csv"
# species <- args[8]
# recons_out_f <- args[9]
decon_out_f <- "/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/st_12_2/Result_Deconv"
# deconv_methoid <- args[11]

st_meta <- read.csv(st_meta, stringsAsFactors=FALSE, row.names=1, check.names=F)
st_data <- read.csv(st_data, row.names=1, check.names=F, stringsAsFactors=FALSE)
sc_data <- read.csv(sc_data, row.names=1, check.names=F, stringsAsFactors=FALSE)
sc_celltype <- read.csv(sc_celltype, check.names=F, stringsAsFactors=FALSE)

label_sc <- sc_celltype$celltype
names(label_sc) <- sc_celltype$cell
label_sc <- as.factor(label_sc)
design.sc <- model.matrix(~label_sc-1)
colnames(design.sc) <- levels(label_sc)
rownames(design.sc) <- rownames(sc_data)
# dge.sc <- deconvSeq::getdge(as.matrix(sc_data), design.sc, ncpm.min=1, nsamp.min=4, method="bin.loess")
# b0.sc <- deconvSeq::getb0.rnaseq(dge.sc, design.sc, ncpm.min=1, nsamp.min=4)
# dge.st <- deconvSeq::getdge(as.matrix(st_data), NULL, ncpm.min=1, nsamp.min=4, method="bin.loess")
# res <- deconvSeq::getx1.rnaseq(NB0=200, b0.sc, dge.st)
# st_coef <- as.matrix(res$x1)

# write.csv(st_coef ,paste0(decon_out_f,'/deconvSeq.csv',sep = "",collapse = NULL),row.names = FALSE)