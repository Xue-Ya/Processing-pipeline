library(SpaTalk)


args = commandArgs()

if (length(args)==0) {
  stop("not enough input", call.=FALSE)
}

count_f <-  args[4]
meta_f <-  args[5]
sc_count <- args[6]
sc_celltype <- args[7]
species <- args[8]
recons_out_f <- args[9]
decon_out_f <- args[10]
deconv_methoid <- args[11]

colData <- read.csv(meta_f, stringsAsFactors=FALSE, row.names=1, check.names=F)
st_counts <- read.csv(count_f, row.names=1, check.names=F, stringsAsFactors=FALSE)
sc_counts <- read.csv(sc_count, row.names=1, check.names=F, stringsAsFactors=FALSE)
sc_celltypes <- read.csv(sc_celltype, row.names=1, check.names=F, stringsAsFactors=FALSE)

obj <- createSpaTalk(st_data = t(as.matrix(st_counts)),
                     st_meta = colData[-4],
                     species = species,
                     if_st_is_sc = F,#
                     spot_max_cell = 10,#
)
                    #  celltype = colData$celltype)

 
# print(ncol(sc_data))
# print(length(sc_celltypes$celltype))

if (deconv_methoid == 2){
  sc_data = (as.matrix(as.integer(unlist(sc_counts))))
}else{
  sc_data = (as.matrix(sc_counts))
}

print(ncol(sc_data))
print(length(sc_celltypes$celltype))

obj <- dec_celltype(obj, method = deconv_methoid, sc_data = sc_data, as.character(sc_celltypes$celltype), use_n_cores = 20,iter_num = 200, if_doParallel = T, dec_result = NULL)

# obj <- find_lr_path(object = obj , lrpairs = lrpairs, pathways = pathways, if_doParallel = T, use_n_cores=20)

# cellname <- unique(colData$celltype)

# for (i in 1:length(cellname)) {
#   try(obj <- dec_cci(object = obj, 
#                      celltype_sender = cellname[i],
#                      celltype_receiver =  cellname[i], 
#                      #  pvalue=0.1, 
#                      if_doParallel = T,  use_n_cores=20))
# }

# n_neis <- c(5, 10, 15, 20, 25, 30)
# pvals <- c(0.01, 0.05, 0.1, 0.15, 0.2)
# min_p <- c(2, 5, 10)
# co_exp <- c(0.05, 0.1, 0.2, 0.5)
# ''
# for (i in 1:length(n_neis)) {
#   for (j in 1:length(pvals)) {
#     for (k in 1:length(min_p)) {
#       for (l in 1:length(co_exp)) {
#         print(paste(pvals[j],n_neis[i], min_p[k], co_exp[l]))
#         try(obj <- dec_cci(object = obj, 
#                     celltype_sender = "Stroma",
#                     celltype_receiver =   "Stroma", 
#                     pvalue=pvals[j], n_neighbor = n_neis[i],
#                     min_pairs = min_p[k], co_exp_ratio = co_exp[l],
#                     if_doParallel = F))
#       }
#     }
#   }
# }

# obj <- dec_cci_all(object = obj, if_doParallel = T,  use_n_cores=20)


# obj <- dec_cci_all(object = obj, pvalue=0.1, if_doParallel = F,  use_n_cores=10)

# celltype_pair <- NULL
# for (i in 1:length(cellname)) {
#     d1 <- data.frame(celltype_sender = rep(cellname[i], length(cellname)), celltype_receiver = cellname,
#         stringsAsFactors = F)
#     celltype_pair <- rbind(celltype_pair, d1)
# }

# tf_path <- NULL
# path_pvalue <- NULL

# for (i in 1:nrow(celltype_pair)) {
#     celltype_sender <- celltype_pair$celltype_sender[i]
#     celltype_receiver <- celltype_pair$celltype_receiver[i]
#     try({obj_lr_path <- get_lr_path(object = obj, celltype_sender = celltype_sender, celltype_receiver = celltype_receiver)
#     tf_path <- rbind(tf_path, obj_lr_path$tf_path)
#     path_pvalue <- rbind(path_pvalue, obj_lr_path$path_pvalue)})
# }

save(obj, file = paste0(recons_out_f,"/spatalk/_spatalk.RData"))
# data_matrix <- as.matrix(obj@data[3],"matrix",sep = "")

# dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/stage2_st_12'
# write.csv(data_matrix ,paste0(out_f,'/spatalk/exp.csv',sep = "",collapse = NULL),row.names = FALSE)
genes <- do.call(rbind,lapply(obj@data[["newdata"]]@Dimnames[1], data.frame))
cells <- do.call(rbind,lapply(obj@data[["newdata"]]@Dimnames[2], data.frame))

write.csv(genes ,paste0(recons_out_f,'/spatalk/genes.csv',sep = "",collapse = NULL),row.names = FALSE)
write.csv(cells ,paste0(recons_out_f,'/spatalk/cells.csv',sep = "",collapse = NULL),row.names = FALSE)

# meta = obj@meta[["newmeta"]]
# write.csv(meta,paste0(out_f,'/spatalk/meta.csv',sep = "",collapse = NULL),row.names = FALSE)

meta = obj@meta[["newmeta"]]
write.csv(meta,paste0(recons_out_f,'/spatalk/meta.csv',sep = "",collapse = NULL),row.names = FALSE)

dense_matrix <- as.matrix(obj@data$newdata)
write.csv(dense_matrix ,paste0(recons_out_f,'/spatalk/exp.csv',sep = "",collapse = NULL),row.names = FALSE)

decon = obj@meta[["rawmeta"]]
write.csv(decon,paste0(decon_out_f ,sep = "",collapse = NULL),row.names = FALSE)

# write.csv(obj@lrpair, paste0(out_f,"/spatalk/_lrpair.csv"), row.names = TRUE)

# write.csv(tf_path, paste0(out_f,"_tfpath.csv"), row.names = TRUE)
# write.csv(path_pvalue, paste0(out_f,"_pathpvalue.csv"), row.names = TRUE)

# load(r_data)


# data_matrix <- as.data.frame(obj@data[["newdata"]]@Dimnames)

# dir <- paste0(dir,'/Result_Repair/spatalk/', sep = "",collapse = NULL)

# write.csv(data_matrix ,paste0(dir,'exp.csv',sep = "",collapse = NULL),row.names = FALSE)




