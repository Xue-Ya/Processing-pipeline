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
out_f <- args[9]

# colData <- read.csv(meta_f, stringsAsFactors=FALSE, row.names=1, check.names=F)
# st_counts <- read.csv(count_f, row.names=1, check.names=F, stringsAsFactors=FALSE)
# sc_counts <- read.csv(sc_count, row.names=1, check.names=F, stringsAsFactors=FALSE)
# sc_celltypes <- read.csv(sc_celltype, row.names=1, check.names=F, stringsAsFactors=FALSE)
obj_path <- paste0(out_f,'/spatalk/_spatalk.RData', sep = "",collapse = NULL)

load(obj_path)
class(obj)
# obj <- createSpaTalk(st_data = t(as.matrix(st_counts)),
#                      st_meta = colData[-4],
#                      species = species,
#                      if_st_is_sc = F,#
#                      spot_max_cell = 10,#
# )
                    #  celltype = colData$celltype)

# obj <- dec_celltype(obj, sc_data = t(as.matrix(sc_counts)), as.character(sc_celltypes$celltype), use_n_cores = 20,iter_num = 500, if_doParallel = T)

obj <- find_lr_path(object = obj , lrpairs = lrpairs, pathways = pathways, if_doParallel = T, use_n_cores=20)

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
''
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

obj <- dec_cci_all(object = obj, if_doParallel = T,  use_n_cores=20)


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

save(obj, file = paste0(out_f,"/spatalk/_spatalk.RData"))
# data_matrix <- as.matrix(obj@data[3],"matrix",sep = "")

# dir = '/home/xuezhengyang/data6/02-deconv_1/Script/FigureData/Dataset1/stage2_st_12'
# write.csv(data_matrix ,paste0(out_f,'/spatalk/exp.csv',sep = "",collapse = NULL),row.names = FALSE)

meta = obj@meta[["newmeta"]]
write.csv(meta,paste0(out_f,'/spatalk/meta.csv',sep = "",collapse = NULL),row.names = FALSE)
write.csv(obj@lrpair, paste0(out_f,"/spatalk/_lrpair.csv"), row.names = TRUE)

# write.csv(tf_path, paste0(out_f,"_tfpath.csv"), row.names = TRUE)
# write.csv(path_pvalue, paste0(out_f,"_pathpvalue.csv"), row.names = TRUE)
