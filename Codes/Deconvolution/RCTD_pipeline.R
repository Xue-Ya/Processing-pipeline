
library(spacexr)
library(Matrix)

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

st_meta <- read.csv(st_meta, stringsAsFactors=FALSE, row.names=1, check.names=F)
st_data <- read.csv(st_data, row.names=1, check.names=F, stringsAsFactors=FALSE)
sc_data <- read.csv(sc_data, row.names=1, check.names=F, stringsAsFactors=FALSE)
sc_celltype <- read.csv(sc_celltype, row.names=1, check.names=F, stringsAsFactors=FALSE)
# convs <- read.csv(conv, row.names=1, check.names=F, stringsAsFactors=FALSE)

cellname <- unique(sc_celltype$celltype)
ref_data <- list()
ref_celltype <- list()
for (i in 1:length(cellname)) {
    ref_celltype1 <- sc_celltype[sc_celltype$celltype == cellname[i], ]
    if (nrow(ref_celltype1) < 25) {
        set.seed(i)
        cell_num <- sample(x = 1:nrow(ref_celltype1), size = 25, replace = T)
        cell_new <- ref_celltype1$cell[cell_num]
        ref_data1 <- sc_data[, cell_new]
        ref_celltype1 <- data.frame(cell = colnames(ref_data1), celltype = cellname[i], stringsAsFactors = F)
    } else{
        ref_data1 <- sc_data[, ref_celltype1$cell]
    }
    ref_data[[i]] <- ref_data1
    ref_celltype[[i]] <- ref_celltype1
}
ref_data1 <- ref_data[[1]]
ref_celltype1 <- ref_celltype[[1]]
for (i in 2:length(cellname)) {
    ref_data1 <- cbind(ref_data1, ref_data[[i]])
    ref_celltype1 <- rbind(ref_celltype1, ref_celltype[[i]])
}
ref_celltype1$cell <- paste0("C", 1:nrow(ref_celltype1))
colnames(ref_data1) <- ref_celltype1$cell
# reference
ref_cell_types <- ref_celltype1$celltype
names(ref_cell_types) <- ref_celltype1$cell
ref_cell_types <- factor(ref_cell_types)
ref_nUMI <- colSums(ref_data1)
ref_refernce <- spacexr::Reference(ref_data1, ref_cell_types, ref_nUMI)
# test data
test_spot_coords <- data.frame(xcoord = st_meta$x, ycoord = st_meta$y)
rownames(test_spot_coords) <- st_meta[,1]
test_spot_nUMI <- colSums(st_data)
test_spot_RCTD <- spacexr::SpatialRNA(test_spot_coords, st_data, test_spot_nUMI)
# run RCTD
myRCTD <- spacexr::create.RCTD(test_spot_RCTD, ref_refernce, max_cores = 10, test_mode = TRUE)
myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = "full")
results <- myRCTD@results
st_coef <- as.matrix(results$weights)

write.csv(st_coef ,paste0(decon_out_f,'/RCTD.csv',sep = "",collapse = NULL),row.names = FALSE)