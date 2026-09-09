## 10x T cells of B6 and M samples; B cells of 3 samples
# QC and total cell clustering
h5read <- Read10X_h5("./filtered_feature_bc_matrix.h5")
raw_data <- CreateSeuratObject(counts = h5read, project = "Seurat", min.cells = 0, min.features = 0)

seurat <- PercentageFeatureSet(raw_data, pattern = MT_name, col.name = "percent.mt")
seurat$orig.ident <- seurat$identy
seurat$identy <- NULL

seurat.data <- subset(seurat, subset =  percent.mt < 10 & nFeature_RNA > 200  & nCount_RNA > 200 )

lm.model <- lm(log10(seurat.data$nFeature_RNA) ~ log10(seurat.data$nCount_RNA))
lm.model$coefficients
seurat.data$valideCells <- log10(seurat.data$nFeature_RNA) > (log10(seurat.data$nCount_RNA) * lm.model$coefficients[2] + 
                                                                (lm.model$coefficients[1] - 0.2))
seurat.data.keep <- seurat.data[,seurat.data$valideCells == "TRUE"]

save_scrublet_data <- function(identities = identities) {
  for (f in identities) {
    seuratdata_subset <- subset(seurat, idents = f)
    
    dir_list <- dir.create(paste("./general_analysis/quality_control/doublet/", f, sep = ""))
    print(dir_list)
    
    counts_matrix <- GetAssayData(seuratdata_subset, assay='RNA', slot='counts')
    writeMM(counts_matrix,  file = paste("./general_analysis/quality_control/doublet/", f, "/count_matrix_sparse.mtx", sep=""))
    
    barcode <- colnames(seuratdata_subset)
    write.table(barcode, file = paste("./general_analysis/quality_control/doublet/",f,"/barcode.tsv", sep=""), 
                quote=F, row.names=F, col.names = F)
  }
} 
Idents(seurat) <- seurat$orig.ident 
identities <- levels(as.factor(seurat$orig.ident))
identities

scrublet <- save_scrublet_data(identities)

read_doublet <- function(identities = identities) {
  for (f in identities) {
    doublet <- read.csv(file=paste("./general_analysis/quality_control/doublet/", f, "/doublet.csv",sep = ""))
    write.table(doublet, file = paste("./general_analysis/quality_control/doublet/total_doublet/",f, "_doublet.txt", sep = ""),
                quote=F, row.names=F)
  }
}
all_doublet <- read_doublet(identities = identities)

{
  a = list.files("./general_analysis/quality_control/doublet/total_doublet") 
  a
  dir = paste("./general_analysis/quality_control/doublet/total_doublet/",a, sep = "")
  n = length(dir)
  n
  merge.doublet.data <- read.table(file = dir[1],header = T, sep = "")
  for (i in 2:n){
    new.doublet.data <- read.table(file = dir[i],header = T, sep = "")
    merge.doublet.data <- rbind(merge.doublet.data,new.doublet.data)
  }
  write.csv(merge.doublet.data, file = "./general_analysis/quality_control/doublet/total_doublet/combined.double.csv", row.names = F)
}

merge.doublet.data <- read.csv("./general_analysis/quality_control/doublet/total_doublet/combined.double.csv")
doublet <- merge.doublet.data$barcode[merge.doublet.data$predicted_doublets == "True"]

seurat.filter <- seurat.data.keep[,!colnames(seurat.data.keep) %in% doublet]

seurat.harmony <- Seurat::NormalizeData(seurat.filter, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = F) %>% 
  RunPCA(npcs = 40, verbose = F)
seurat.harmony <- seurat.harmony %>% 
  RunHarmony("orig.ident", plot_convergence = T, reduction.save = "harmony")

seurat_harmony <- RunUMAP(seurat.harmony, n.components = 2, reduction = "harmony",
                          min.dist = 0.3, dims = 1:40, seed.use = 3)
seurat_harmony <- FindNeighbors(seurat_harmony, reduction = "harmony", dims = 1:40)
seurat.cluster <- FindClusters(seurat_harmony,resolution=0.3)

# T cells of sample M
M.cluster <- Seurat::NormalizeData(M, verbose = FALSE)
M.cluster <- CellCycleScoring(M.cluster, s.features = m.s.genes, 
                              g2m.features = m.g2m.genes, set.ident = TRUE)
M.cluster <- FindVariableFeatures(M.cluster, selection.method = "vst", nfeatures = 2000) 

M.cluster <- ScaleData(M.cluster, vars.to.regress = c("S.Score", "G2M.Score"), 
                       features = rownames(M.cluster))

M.clustering <- RunPCA(M.cluster, npcs = 30, verbose = F)

M.clustering <- RunUMAP(M.clustering, n.components = 2, min.dist = 0.3, dims = 1:30, seed.use = 3)
M.clustering <- FindNeighbors(M.clustering, dims = 1:30)
M.clustering <- FindClusters(M.clustering)

# T cells of sample B6
B6.reclustering <- Seurat::NormalizeData(B6, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

B6.reclustering <- ScaleData(B6.reclustering, verbose = F) %>% 
  RunPCA(npcs = 20, verbose = F)

B6.reclustering <- RunUMAP(B6.reclustering, n.components = 2, min.dist = 0.3, dims = 1:20, seed.use = 3)
B6.reclustering <- FindNeighbors(B6.reclustering, dims = 1:20)
B6.reclustering <- FindClusters(B6.reclustering)


# B cells of 3 samples
Bcells.reclustering <- Seurat::NormalizeData(Bcells, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

Bcells.reclustering <- ScaleData(Bcells.reclustering, verbose = F) %>% 
  RunPCA(npcs = 15, verbose = F)

system.time({Bcells.reclustering <- RunHarmony(Bcells.reclustering, group.by.vars = "orig.ident")})
Bcells.reclustering <- RunUMAP(Bcells.reclustering, reduction = "harmony", min.dist = 0.3, dims = 1:15)
Bcells.reclustering <- FindNeighbors(Bcells.reclustering, dims = 1:15)
Bcells.reclustering <- FindClusters(Bcells.reclustering, resolution = 0.8)

MGs <- FindAllMarkers(Bcells.reclustering, verbose = FALSE, min.cells.group = 5, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)










