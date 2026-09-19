# QC
library(Seurat)
library(DoubletFinder)

work_dir="./02.matrix_in"
setwd(work_dir)

sample.mex.list <- list.files(path=work_dir)
sample.mex.list

dir <- paste(work_dir,sample.mex.list,sep="/")
dir

rawdata.list <- list()

for(i in dir){
  filename <- list.files(path=i)
  filepath <- paste(i,filename,sep="/")
  raw.matrix <- ReadMtx(mtx = filepath[3],
                        cells = filepath[1],
                        features = filepath[2])
  rawdata.list[[i]] <- raw.matrix
}

names(rawdata.list) <- sample.mex.list
lapply(rawdata.list,dim)

obj.list <- lapply(X=names(rawdata.list),FUN = function(x){
  x <- CreateSeuratObject(counts = rawdata.list[[x]],
                          project = x,
                          min.cells = 3,
                          min.features = 200)
})

names(obj.list) <- names(rawdata.list)
obj.list

for (i in names(obj.list)){
  obj.list[[i]][["percent.mt"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^mt-")
  obj.list[[i]][["orig.ident"]] <- i
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], selection.method = "vst", nfeatures = 2000)
  obj.list[[i]] <- ScaleData(obj.list[[i]])
  obj.list[[i]] <- RunPCA(obj.list[[i]])
  obj.list[[i]] <- FindNeighbors(obj.list[[i]], dims = 1:10)
  obj.list[[i]] <- FindClusters(obj.list[[i]], resolution = 0.5)
  obj.list[[i]] <- RunUMAP(obj.list[[i]], dims = 1:10)
  
  sweep.res.list_kidney <- paramSweep(obj.list[[i]], PCs = 1:10, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  
  pK_bcmvn <- as.numeric(as.character(bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]))
  annotations <- obj.list[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)     
  
  nExp_poi <- round(0.055*nrow(obj.list[[i]]@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  obj.list[[i]] <- doubletFinder(obj.list[[i]], PCs = 1:10, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  obj.list[[i]] <- doubletFinder(obj.list[[i]], PCs = 1:10, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  print(colnames(obj.list[[i]]@meta.data))
  colnames(obj.list[[i]]@meta.data) <- c('orig.ident','nCount_RNA','nFeature_RNA','percent.mt','RNA_snn_res.0.5','seurat_clusters','pANN','DF','pANN_adj','DF_adj')
}

head(obj.list$`after_LN`@meta.data)

# 多样本合并
merge.obj <- merge(x = obj.list[[1]],
                   y = obj.list[2:length(obj.list)],
                   add.cell.ids = c("a-1","a-2","a-3","a-4")
)
merge.obj

head(merge.obj@meta.data)

table(merge.obj$DF)
table(merge.obj$DF_adj)

library(ggplot2)
library(dplyr)

work_dir <- "./03.out"
setwd(work_dir)
setwd("1.QC")

# Violin plot
options(repr.plot.width=12, repr.plot.height=6)
qc_feature <- c("nFeature_RNA","nCount_RNA","percent.mt")
p2 <- VlnPlot(merge.obj,features = qc_feature, ncol=3,pt.size = 0,group.by = 'orig.ident')
p2
ggsave(p2,file="2_VlnPlot_BeforeQC.pdf",width = 12,height=6)

p3 <- FeatureScatter(merge.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")
p4 <- FeatureScatter(merge.obj, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")
p3 + p4 
ggsave(p3+p4,file="3_FeatureScatter_BeforeQC.pdf",width = 10,height=5)

merge.obj$log10GenesPerUMI <- log10(merge.obj$nFeature_RNA) / log10(merge.obj$nCount_RNA)
options(repr.plot.width=6, repr.plot.height=7)
p8 <- merge.obj@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  facet_wrap("orig.ident ~.",ncol=1,strip.position='left') +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.88)
p8
ggsave(p8,file="4_density_log10GenesPerUMI_BeforeQC.pdf",width = 6,height=7)

merge.obj <- NormalizeData(merge.obj, verbose = TRUE)
merge.obj <- FindVariableFeatures(merge.obj, selection.method = "vst", nfeatures = 2000)
merge.obj <- ScaleData(merge.obj) 
merge.obj <- RunPCA(merge.obj, features = VariableFeatures(object = merge.obj))
merge.obj <- FindNeighbors(merge.obj, dims=1:30)
merge.obj <- FindClusters(merge.obj, resolution = 0.5, cluster.name = "cluster.unintegrated")
merge.obj <- RunUMAP(merge.obj, dims = 1:30, reduction.name = "umap.unintergrated")

merge.obj@meta.data[,"DF_hi.lo"] <- merge.obj@meta.data$DF
merge.obj@meta.data$DF_hi.lo[which(merge.obj@meta.data$DF_hi.lo == "Doublet" & merge.obj@meta.data$DF_adj == "Singlet")] <- "Doublet-Low Confidience"
merge.obj@meta.data$DF_hi.lo[which(merge.obj@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"

merge.obj <- subset(merge.obj,subset =  nCount_RNA>300 & nFeature_RNA>200 & nFeature_RNA <5000 & percent.mt < 15 &log10GenesPerUMI > 0.88)
merge.obj <- subset(merge.obj, subset= DF_hi.lo %in% c("Singlet"))

# scRNA-seq cluster
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)

# 读入QC后的细胞
bc <- read.csv("/storage/dongchenLab/linmeng/41_yuxiangyu/m1/03.out/1.QC/8_meta.csv",row.names = 'X')

k=0
for (i in names(obj.list)){
  k <- k+1
  obj.list[[i]][["percent.mt"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^mt-") 
  obj.list[[i]][["orig.ident"]] <- i
  j=paste0("a-",k)
  obj.list[[i]] <- RenameCells(obj.list[[i]],add.cell.id = j)
  obj.list[[i]] <- subset(obj.list[[i]], cells = rownames(bc)) 
}

# 多样本合并
merge.obj <- merge(x = obj.list[[1]],
                   y = obj.list[2:length(obj.list)]
)

merge.obj <- NormalizeData(merge.obj)
merge.obj <- FindVariableFeatures(merge.obj,selection.method = "vst", nfeatures = 2000)

# 降维聚类
VariableFeatures(merge.obj) <- VariableFeatures(merge.obj)[!grepl("^Tr[a/b/g/d][v/d/j/c]",VariableFeatures(merge.obj))] 

library(homologene)
ts <- human2mouse(cc.genes.updated.2019$s.genes)
s.genes <- ts$mouseGene
ts <- human2mouse(cc.genes.updated.2019$g2m.genes)
g2m.genes <- ts$mouseGene

merge.obj<-JoinLayers(merge.obj)

merge.obj <- CellCycleScoring(merge.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) 

system.time({
  all.genes <- rownames(merge.obj)
  merge.obj <- ScaleData(merge.obj, features = all.genes,vars.to.regress = c("S.Score", "G2M.Score")) 
})

merge.obj <- RunPCA(merge.obj,features = VariableFeatures(object = merge.obj))
merge.obj <- RunHarmony(merge.obj,group.by.vars = "orig.ident",plot_convergence = TRUE)
merge.obj <- FindNeighbors(merge.obj, dims = 1:20,reduction = "harmony")
merge.obj <- FindClusters(merge.obj, resolution = 0.4)
merge.obj <- RunUMAP(merge.obj, dims = 1:20,reduction = "harmony")

markers <- FindAllMarkers(merge.obj, only.pos = TRUE)
write.csv(markers,file="1.marker_cluster.csv")
markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  top_n(n=10,wt=avg_log2FC) -> top_mk


# Cluster subset&recluster
sub <- subset(df,seurat_clusters %in% "0") #or C3
sub <- NormalizeData(sub)
sub <- FindVariableFeatures(sub,selection.method = "vst", nfeatures = 2000)
VariableFeatures(sub) <- VariableFeatures(sub)[!grepl("^Tr[a/b/g/d][v/d/j/c]",VariableFeatures(sub))]
library(homologene)
ts <- human2mouse(cc.genes.updated.2019$s.genes)
s.genes <- ts$mouseGene
ts <- human2mouse(cc.genes.updated.2019$g2m.genes)
g2m.genes <- ts$mouseGene
sub <- CellCycleScoring(sub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

all.genes <- rownames(sub)
sub <- ScaleData(sub, vars.to.regress = c("S.Score", "G2M.Score")) 
sub <- RunPCA(sub,features = VariableFeatures(object = sub))
options(repr.plot.width=5, repr.plot.height=5)
ElbowPlot(sub,ndims =50)

sub <- RunHarmony(sub,group.by.vars = "orig.ident",plot_convergence = TRUE)
sub <- FindNeighbors(sub, dims = 1:10,reduction = "harmony")
sub <- FindClusters(sub, resolution = 0.3)
sub <- RunUMAP(sub, dims = 1:10,reduction = "harmony")

options(repr.plot.width=5, repr.plot.height=5)
p<-DimPlot(sub, reduction = "umap", pt.size = 0.1,group.by = "seurat_clusters",raster=FALSE,label = T)
p
ggsave(p,file="DimPlot_cluster.pdf")

# Celltype anno
meta <- df@meta.data
meta$celltype <- "no"

meta2 <- read.csv("/4.sub_C0/1.meta.csv",row.names = 'X')
meta[rownames(meta2),]$celltype <- meta2$celltype
meta[meta$seurat_clusters %in% c('1','2','10'),]$celltype <- "Naïve CD4+ T cell"
meta2 <- read.csv("/4.sub_C3/3.meta_final.csv",row.names = 'X')
meta[rownames(meta2),]$celltype <- meta2$celltype

meta[meta$seurat_clusters %in% c('4'),]$celltype <- "Tfh"
meta[meta$seurat_clusters %in% c('5'),]$celltype <- "IFN induced CD4+ T cell"
meta[meta$seurat_clusters %in% c('6'),]$celltype <- "Early priming CD4+ T cell"
meta[meta$seurat_clusters %in% c('7'),]$celltype <- "Mbd4+ CD4+ T cell"
meta[meta$seurat_clusters %in% c('8'),]$celltype <- "Cytotoxic CD4+ T cell"
meta[meta$seurat_clusters %in% c('9'),]$celltype <- "Proliferating CD4+ T cell"

df@meta.data <- meta

options(repr.plot.width=8.2,repr.plot.height=5.7)
p<-DimPlot(df,group.by = 'celltype',cols = col_celltype, raster=FALSE)
p
ggsave(p,file="3.DimPlot_celltype.pdf",width=8.2,height=6)

options(repr.plot.width=14.2,repr.plot.height=4)
p<-DimPlot(df,group.by = 'celltype',cols = col_celltype, split.by = 'orig.ident',raster=FALSE)
p
ggsave(p,file="3.DimPlot_celltype_split.pdf",width=14.2,height=4)

sub <-subset(df,orig.ident %in% c('MN-1','MN-2'))
p<-DimPlot(sub,group.by = 'celltype',cols = col_celltype, raster=FALSE)
p
ggsave(p,file="3.DimPlot_celltype_MN.pdf",width=8.5,height=6)

genes_plot = c('Cxcr5','Tox2','Bcl6')
p <- FeaturePlot(sub,features = genes_plot,order = TRUE,raster=FALSE,ncol = 3)
ggsave(p,file="3.FeaturePlot_MN_Fig5B.pdf",width=16,height=5)

genes_plot = c('Ifng','Il17a','Tbx21','Rorc')
p <- FeaturePlot(sub,features = genes_plot,order = TRUE,raster=FALSE,ncol = 2)
ggsave(p,file="3.FeaturePlot_MN_SuppFig5A.pdf",width=10,height=9)

# cell proportion
Idents(df)<-factor(df$celltype,levels = names(table(df$celltype)))
options(repr.plot.width=8, repr.plot.height=6)
cell.prop<-as.data.frame(prop.table(table(Idents(df), df$orig.ident)))
colnames(cell.prop)<-c("celltype","sample","proportion")

p <- ggplot(cell.prop,aes(sample,proportion,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values = col_celltype)+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))
p
ggsave(p,file="4.cell_proportion_sample_celltype.pdf",width = 8,height = 6)

# dotplot
gene_list <- list(IFN,Naive,Stem,Effector_T,Tfh,Early_p,Th17,Effector_M,Mbd4,Th1,Prolif,Cytotoxic,HSP)
names(gene_list) <- c("IFN induced\nCD4+ T cell","Naïve\nCD4+ T cell","Stem-like Treg","Effector Treg","Tfh",
                      "Early priming\nCD4+ T cell","Th17","Effector Memory\nCD4+ T cell","Mbd4+\nCD4+ T cell",
                      "Th1","Proliferating\nCD4+ T cell","Cytotoxic\nCD4+ T cell","HSP high\nCD4+ T cell")

df <- SetIdent(df, value = 'celltype')

options(repr.plot.width=25, repr.plot.height=7)
p<-DotPlot(object = df, features=gene_list) + 
  scale_size(range=c(0,12))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=15),
        axis.text.y = element_text(size=15), # 横坐标标签倾斜
        strip.text.x = element_text(colour = "black", face = "bold",size=15),
        strip.background = element_rect(colour = "white", fill = "lightgrey"))
p

pdf(file="1.Dotplot_all.pdf",width=60,height=7)
grid.newpage()
dev.off()


# TCR mapping
library(scRepertoire)
library(Seurat)
library(dplyr)
library(stringr)

contig.list <- loadContigs(input = "./02.TCR",format = "BD")
airr <- contig.list[[1]]
airr <- subset(airr,chain %in% c("TRA","TRB")) # 筛选chain为TRA和TRB的
airr <- airr[str_detect(airr$v_gene,"TR[A,B]") | airr$v_gene=="",] # v_gene只是TRA和TRB的，去掉掺杂的TRD，TRG or BCR相关基因；
airr <- airr[str_detect(airr$d_gene,"TR[A,B]") | airr$d_gene=="",] # d_gene只是TRA和TRB的，去掉掺杂的TRD，TRG or BCR相关基因；
airr <- airr[str_detect(airr$j_gene,"TR[A,B]") | airr$j_gene=="",] # j_gene只是TRA和TRB的，去掉掺杂的TRD，TRG or BCR相关基因；
airr <- airr[str_detect(airr$c_gene,"Tr[a,b]") | airr$c_gene=="",] # c_gene只是TRA和TRB的，去掉掺杂的TRD，TRG or BCR相关基因；
contig.list[[1]] <- airr

combined.TCR <- combineTCR(contig.list, 
                           samples = c("a-1","a-2","a-3","a-4"),
                           removeNA = FALSE,
                           removeMulti = FALSE,
                           filterMulti = TRUE, 
                           filterNonproductive = TRUE) 

df <- merge.obj
df <- combineExpression(combined.TCR, 
                        df, 
                        cloneCall="aa", 
                        chain = "both", 
                        proportion = FALSE, 
                        cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))


options(repr.plot.width=25, repr.plot.height=5)
DimPlot(df, group.by = "cloneSize",split.by = "cloneSize",order=TRUE,raster=FALSE) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))

# TCR expansion
meta <- df@meta.data
meta <- meta[!is.na(meta$CTaa),]
meta <- meta[!grepl("^NA_|_NA$",meta$CTaa),]

meta<-subset(meta,orig.ident %in% c('LN-1','LN-2')) # or 'MN-1','MN-2'

cell.prop <- as.data.frame(prop.table(table(meta$cloneSize,meta$celltype)))
colnames(cell.prop)<-c("cloneSize","celltype","proportion")
cell.prop$cloneSize <- factor(cell.prop$cloneSize, levels=rev(c('Hyperexpanded (100 < X <= 642)','Large (20 < X <= 100)',
                                                                'Medium (5 < X <= 20)','Small (1 < X <= 5)','Single (0 < X <= 1)')))

options(repr.plot.width=7, repr.plot.height=5)
p <- ggplot(cell.prop,aes(celltype,proportion,fill=cloneSize))+
  geom_bar(stat="identity",position="fill")+
  scale_y_continuous(labels = scales::percent_format(accuracy=1)) +    
  ggtitle("LN")+
  xlab("")+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title="Clonal Group"))+
  scale_fill_manual(values = col_cloneSize)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size=15), 
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(angle = 70,vjust = 1,hjust = 1,size=13) 
  )
p
ggsave(p,file="3.TCR_expansion_LN.pdf",width=7,height=5)


# LN、MN Tfh TCR analysis
library(scRepertoire)
library(circlize)
library(scales)
library(stringr)

library(ggalluvial)
library(ggplot2)
library(grid)
library(ggbreak)

meta <- meta[meta$celltype %in% 'Tfh',]
meta$group <- "LN"
meta[meta$orig.ident %in% c("MN-1","MN-2"),]$group <- "MN"
meta$barcode <- rownames(meta)
meta <- meta[!is.na(meta$CTaa),]
meta <- meta[!grepl("^NA_|_NA$",meta$CTaa),] 

meta1 <- subset(meta,group %in% c("LN")) 
meta2 <- subset(meta,group %in% c("MN")) 

prop_eff <- table(meta1$CTaa)
prop_eff <- as.data.frame(prop_eff)
prop_eff <- prop_eff[order(prop_eff$Freq,decreasing = TRUE),]
rownames(prop_eff) <- NULL
colnames(prop_eff) <- c('clonetype','Freq')
prop_eff$celltype2 <- "LN"

prop_stm <- table(meta2$CTaa)
prop_stm <- as.data.frame(prop_stm)
prop_stm <- prop_stm[order(prop_stm$Freq,decreasing = TRUE),]
rownames(prop_stm) <- NULL
colnames(prop_stm) <- c('clonetype','Freq')
prop_stm$celltype2 <- "MN"

prop_eff$pro <- prop.table(prop_eff$Freq) 
prop_stm$pro <- prop.table(prop_stm$Freq)  

prop_eff$freq_abs <- prop_eff$Freq/4831 
prop_stm$freq_abs <- prop_stm$Freq/1962 

data <- rbind(prop_eff,prop_stm)
data$celltype2 <- factor(data$celltype2,levels=c('LN','MN'))

overlap <- merge(prop_eff,prop_stm,by = 'clonetype')
data1 <- data[data$clonetype %in% overlap$clonetype,]
data1$id <- 'yes'
data2 <- data[!(data$clonetype %in% overlap$clonetype),]
data2$id <- 'no'
data <- rbind(data1,data2)
data$id <- factor(data$id,levels=c('yes','no'))

options(repr.plot.width=6, repr.plot.height=6)
p<-ggplot(data,aes(x = celltype2, y = pro, alluvium = clonetype, stratum = id))+
  scale_y_break(c(0.03, 0.94))+ 
  geom_alluvium(width = 0.3, aes(fill = clonetype)) +
  geom_stratum(width = 0.3, 
               fill=c('#54b4b1','#54b4b1','#a96fa9','#a96fa9'),
               color=c('#54b4b1','#54b4b1','#a96fa9','#a96fa9'))+ 
  coord_cartesian(ylim = c(0.9, 1))+
  ylab("proportion(%)") +
  labs(title="Tfh")+
  theme_classic() +
  theme(plot.title=element_text(face="bold",size=24,hjust=0.5,vjust=0.5),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size =18),
        axis.text.y = element_text(size =18),
        axis.title.y = element_text(size =20),
        legend.position = "none" 
  )
p
ggsave(p,file='2.1_TCRsharing_Tfh_LN_MN_1.pdf',width=6, height=6,onefile = FALSE)

# TCR Diversity
library(immunarch)

meta <- meta[!is.na(meta$CTaa),]
meta <- meta[!grepl("^NA_|_NA$",meta$CTaa),]

md_table <- data.frame(
  filename = c("LN-1_VDJ_Dominant_Contigs_AIRR.tsv",
               "LN-2_VDJ_Dominant_Contigs_AIRR.tsv",
               "MN-1_VDJ_Dominant_Contigs_AIRR.tsv",
               "MN-2_VDJ_Dominant_Contigs_AIRR.tsv"),
  sample   = c("LN-1", "LN-2", "MN-1","MN-2"),
  group    = c("LN", "LN", "MN","MN")
)

schema <- make_receptor_schema(features = c("cdr3_aa", "v_call","d_call","j_call","c_call"), chains = c("TRA", "TRB"))

idata <- read_repertoires(path = "/02.TCR_rename/*tsv", 
                          schema = schema, 
                          metadata = md_table,
                          barcode_col = "cell_id",
                          locus_col = "locus",
                          umi_col = "umi_count", 
                          preprocess = make_default_preprocessing("airr"), 
                          repertoire_schema = NULL)

idata <- annotate_barcodes(idata, tibble(meta), "barcode")

idata <- agg_repertoires(idata, c('group.y','orig.ident'))

idata <- filter_barcodes(idata, meta$barcode, keep_repertoires = TRUE)

data <- airr_diversity_dxx(
  idata,
  perc = 50,
  autojoin = getOption("immundata.autojoin", TRUE),
  format = c("long", "wide")
)

data$group <- data$group.y

p<-ggplot(data,aes(x=group.y,y=dxx))+
  geom_boxplot(aes(color=group),linewidth=1)+
  scale_color_manual(values = c("LN" = "#54b4b1", "MN" = "#a96fa9"))+
  geom_jitter(color="#696969",size=3,width=0.05)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 18,face='bold'),
        axis.text.y = element_text(size = 18,face='bold'),
        axis.line.x = element_line(color = "black", size = 1, linetype = "solid"),
        axis.line.y = element_line(color = "black", size = 1, linetype = "solid"),
        plot.title = element_text(size = 20,face='bold')
  ) +
  ggtitle("D50")
p
ggsave(p,file="2.2_D50.pdf",width=6,height=6)

# Th1 or Th17 TCR sharing
meta <- meta[meta$celltype %in% 'Th1',] # or Th17
meta$group <- "LN"
meta[meta$orig.ident %in% c("MN-1","MN-2"),]$group <- "MN"
meta$barcode <- rownames(meta)

meta <- meta[!is.na(meta$CTaa),]
meta <- meta[!grepl("^NA_|_NA$",meta$CTaa),] 

meta1 <- subset(meta,group %in% c("LN")) 
meta2 <- subset(meta,group %in% c("MN"))

prop_eff <- table(meta1$CTaa) 
prop_eff <- as.data.frame(prop_eff)
prop_eff <- prop_eff[order(prop_eff$Freq,decreasing = TRUE),]
rownames(prop_eff) <- NULL
colnames(prop_eff) <- c('clonetype','Freq')
prop_eff$celltype2 <- "LN"

prop_stm <- table(meta2$CTaa)
prop_stm <- as.data.frame(prop_stm)
prop_stm <- prop_stm[order(prop_stm$Freq,decreasing = TRUE),]
rownames(prop_stm) <- NULL
colnames(prop_stm) <- c('clonetype','Freq')
prop_stm$celltype2 <- "MN"

prop_eff$pro <- prop.table(prop_eff$Freq) 
prop_stm$pro <- prop.table(prop_stm$Freq)

prop_eff$freq_abs <- prop_eff$Freq/568  # Th17 1038
prop_stm$freq_abs <- prop_stm$Freq/10971 # Th17 8928

data <- rbind(prop_eff,prop_stm)
data$celltype2 <- factor(data$celltype2,levels=c('LN','MN'))

overlap <- merge(prop_eff,prop_stm,by = 'clonetype')
data1 <- data[data$clonetype %in% overlap$clonetype,]
data1$id <- 'yes'
data2 <- data[!(data$clonetype %in% overlap$clonetype),]
data2$id <- 'no'
data <- rbind(data1,data2)
data$id <- factor(data$id,levels=c('yes','no'))

options(repr.plot.width=6, repr.plot.height=6)
p<-ggplot(data,aes(x = celltype2, y = pro, alluvium = clonetype, stratum = id))+
  scale_y_break(c(0.1, 0.75))+ 
  geom_alluvium(width = 0.3, aes(fill = clonetype)) +
  geom_stratum(width = 0.3, 
               fill=c('#54b4b1','#54b4b1','#a96fa9','#a96fa9'),
               color=c('#54b4b1','#54b4b1','#a96fa9','#a96fa9'))+
  coord_cartesian(ylim = c(0.75, 1))+ 
  ylab("proportion(%)") +
  labs(title="Th1")+
  theme_classic() +
  theme(plot.title=element_text(face="bold",size=24,hjust=0.5,vjust=0.5),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size =18),
        axis.text.y = element_text(size =18),
        axis.title.y = element_text(size =20),
        legend.position = "none" 
  )
p
ggsave(p,file='2.1_TCRsharing_Th1_LN_MN_1.pdf',width=6, height=6,onefile = FALSE)

# B cells and Tfh cell-cell interaction
library(Seurat)
library(ggplot2)
library(grid)
library(viridis)
library(CellChat)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggraph)

load("Seuratobject_B_cells.RData")
SeuObj_B$celltype <- "Total B cell"

seurat_object

# sub: MN （sub: LN同理）
j=c("MN-1","MN-2")

for(i in c("35-55","B6","rhu-1")){
  seurat_object=cp
  seurat_object <- subset(seurat_object, orig.ident %in% c(i,j))
  seurat_object$celltype <- droplevels(seurat_object$celltype) 
  name <- paste0("MN_",i)
  setwd(work_dir)
  dir.create(name)
  setwd(name) 
  # cellchat analysis
  cellchat <- createCellChat(object = seurat_object, group.by = "celltype", assay = "RNA")
  CellChatDB <- CellChatDB.mouse
  CellChatDB.use <- subsetDB(CellChatDB,
                             search=c("CXCL10_CXCR3","CXCL16_CXCR6",
                                      "CD86_CD28","CD86_CTLA4","ICOSL_ICOS",
                                      "ICOSL_CD28","ICOSL_CTLA4","LGALS9_CD45",
                                      "CXCL13_CXCR5"),
                             key="interaction_name",
                             non_protein = FALSE)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat,min.cells = 10,thresh.p =0.05)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, 
                                type = "truncatedMean", 
                                trim=0)
  df.net <- subsetCommunication(cellchat,thresh=NULL) 
  write.csv(df.net,file="2.1.cellchat_net.csv")
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  col_celltype = scPalette(length(levels(cellchat@idents)))
  names(col_celltype) <- levels(cellchat@idents)
  
  pdf("2.1.netVisual_circle_count.pdf",width=7,height=7)
  p <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                        weight.scale = T, 
                        arrow.size = 0.5,
                        title.name = "Number of interactions",
                        vertex.label.cex  = 0.01,
                        vertex.label.color = "white",
                        margin = 0.5,
                        color.use = col_celltype) 
  legend(x=-1.95, y=-1.3, names(col_celltype),
         pch=16, 
         col=col_celltype, 
         pt.cex=0.8, 
         cex=0.52, 
         bty="n",ncol=5)
  dev.off()
  
  pdf("2.1.netVisual_circle_weight.pdf",width=7,height=7)
  p <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                        weight.scale = T, 
                        arrow.size = 0.5,
                        title.name = "Interaction weights/strength",
                        vertex.label.cex  = 0.01,
                        vertex.label.color = "white",
                        margin = 0.5,
                        color.use = col_celltype)
  legend(x=-1.95, y=-1.3, names(col_celltype),
         pch=16, 
         col=col_celltype, 
         pt.cex=0.8,  
         cex=0.52, 
         bty="n",ncol=5)
  dev.off()
  
  # counts: source
  mat <- cellchat@net$count
  col_celltype <- scPalette(nrow(mat))
  names(col_celltype) <- rownames(mat)
  
  pdf("2.1.cellchat_net_count_source.pdf",width = 20, height = 10)
  par(mfrow = c(2,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ] # source
    netVisual_circle(mat2, vertex.weight = groupSize, 
                     weight.scale = T, 
                     arrow.size = 0.8,
                     title.name = rownames(mat)[i],
                     vertex.label.cex  = 0.01,
                     vertex.label.color = "white",
                     color.use = col_celltype)
    legend(x=-1.5, y=-1.1, names(col_celltype),
           pch=16, 
           col=col_celltype, 
           pt.cex=1,  
           cex=0.6, 
           bty="n",ncol=5)
  }
  dev.off()
  
  # counts: target
  mat <- cellchat@net$count
  pdf("2.1.cellchat_net_count_target.pdf",width = 20, height = 10)
  par(mfrow = c(2,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[,i] <- mat[,i] # target
    netVisual_circle(mat2, vertex.weight = groupSize, 
                     weight.scale = T, 
                     arrow.size = 0.8,
                     title.name = rownames(mat)[i],
                     vertex.label.cex  = 0.01,
                     vertex.label.color = "white",
                     color.use = col_celltype)
    legend(x=-1.5, y=-1.1, names(col_celltype),
           pch=16, 
           col=col_celltype, 
           pt.cex=1,  
           cex=0.6, 
           bty="n",ncol=5)
  }
  dev.off()
  
  # weight: source
  mat <- cellchat@net$weight
  col_celltype <- scPalette(nrow(mat))
  names(col_celltype) <- rownames(mat)
  
  pdf("2.1.cellchat_net_weight_source.pdf",width = 20, height = 10)
  par(mfrow = c(2,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ] # source
    netVisual_circle(mat2, vertex.weight = groupSize, 
                     weight.scale = T, 
                     arrow.size = 0.8,
                     title.name = rownames(mat)[i],
                     vertex.label.cex  = 0.01,
                     vertex.label.color = "white",
                     color.use = col_celltype)
    legend(x=-1.5, y=-1.1, names(col_celltype),
           pch=16, 
           col=col_celltype, 
           pt.cex=1,  
           cex=0.6, 
           bty="n",ncol=5)
  }
  dev.off()
  
  # weight: target
  mat <- cellchat@net$weight
  pdf("2.1.cellchat_net_weight_target.pdf",width = 20, height = 10)
  par(mfrow = c(2,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[,i] <- mat[,i] # target
    netVisual_circle(mat2, vertex.weight = groupSize, 
                     weight.scale = T, 
                     arrow.size = 0.8,
                     title.name = rownames(mat)[i],
                     vertex.label.cex  = 0.01,
                     vertex.label.color = "white",
                     color.use = col_celltype)
    legend(x=-1.5, y=-1.1, names(col_celltype),
           pch=16, 
           col=col_celltype, 
           pt.cex=1,  
           cex=0.6, 
           bty="n",ncol=5)
  }
  dev.off()
  
  p<-netVisual_bubble(cellchat, remove.isolate = FALSE)
  ggsave(p,file="2.1.netVisual_bubble.pdf",width=49,height=15)
}

# Interaction dotplot
# ct_list[4:6]: MN vs 3 B samples
data <- c()
for(j in ct_list[4:6]){
  cellchat <- readRDS(j)
  dt <- subsetCommunication(cellchat, thresh = 1.1) 
  dt <- subset(dt, source %in% c("Total B cell",'Tfh'))
  dt <- subset(dt, target %in% c("Total B cell",'Tfh'))
  dtp <- c()
  for(i in 1:nrow(sig)){
    dtp <- rbind(dtp,dt[dt$ligand==sig[i,1] & dt$receptor==sig[i,2],])
    dtp <- rbind(dtp,dt[dt$ligand==sig[i,2] & dt$receptor==sig[i,1],])
  }
  group_name <- strsplit(j,"/")[[1]][1]
  dtp$xname <- paste0(group_name,": ",dtp$source," -> ",dtp$target)
  dtp <- dtp[dtp$source!=dtp$target,] 
  data <- rbind(data,dtp)
}
data$pval_sig =1
data[data$pval>0.05,]$pval_sig <- '1'
data[data$pval<=0.05 & data$pval>0.01,]$pval_sig = '2' 
data[data$pval<=0.01,]$pval_sig = '3'
data$pval_sig <- as.numeric(data$pval_sig)
values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p <= 0.05", "p <= 0.01")

data2 <- subset(data, source %in% "Total B cell")
data2$pval_sig =1
data2[data2$pval>0.05,]$pval_sig <- '1' 
data2[data2$pval<=0.05 & data2$pval>0.01,]$pval_sig = '2' 
data2[data2$pval<=0.01,]$pval_sig = '3'
data2$pval_sig <- as.numeric(data2$pval_sig)

data2$interaction_name_2 <- factor(data2$interaction_name_2,
                                   levels = rev(c('Cd86  - Cd28','Cd86  - Ctla4','Icosl  - Icos',
                                                  'Cxcl16  - Cxcr6','Cxcl10  - Cxcr3'))) 
p<-ggplot(data2, aes(x = xname, y = interaction_name_2, 
                     color = prob, size = pval_sig)) +
  geom_point(pch = 16) + 
  theme_linedraw() + theme(panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                   vjust = 0.5), axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "bottom")+
  scale_radius(range = c(min(data2$pval_sig), max(data2$pval_sig)), 
               breaks = sort(unique(data2$pval_sig)), labels = names(values)[values %in% sort(unique(data2$pval_sig))], name = "p-value")+
  scale_colour_gradientn(colors = colorRampPalette(c('#5e4fa2','#3288bd','#66c2a5','#abdda4','#e6f598','#fee08b','#fdae61','#f46d43','#d53e4f','#9e0142'))(99), 
                         na.value = "white", limits = c(quantile(data2$prob, 
                                                                 0, na.rm = T), quantile(data2$prob, 1, na.rm = T)), 
                         breaks = c(quantile(data2$prob, 0, na.rm = T), quantile(data2$prob, 
                                                                                 1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5, 
                                                                                                                                                            title = "Commun. Prob."))+
  theme(text = element_text(size = 11), plot.title = element_text(size = 10)) + 
  theme(legend.title = element_text(size = 9), legend.text = element_text(size = 8))
p
ggsave(p,file="1.dotplot_MN_Bcell_source_new.pdf",width=5,height=4.5)

# B cells Dotplot
gene_list <- c('Cd86','Icosl','Cxcl16','Cxcl10','Cxcl13')

options(repr.plot.width=7,repr.plot.height=3.5)
p<-DotPlot(object = SeuObj_B, features=gene_list,group.by = 'orig.ident',dot.scale =8,scale = TRUE)
p
ggsave(p,file="2.DotPlot_B_split_by_sample.pdf",width=7,height=3.5)

# B cells Vlnplots
B_gene <- c("Icosl","Cxcl13","Cd86")

for(i in B_gene){
  plot_data <- FetchData(SeuObj_B, vars = c(get('i'),"celltype","orig.ident"),layer = 'data')
  colnames(plot_data) <- c("gene","celltype","orig.ident")
  p<-ggplot(plot_data, aes(x = celltype, y = gene, fill = orig.ident)) +
    geom_violin(scale = "width", trim = TRUE) + 
    facet_grid(orig.ident ~ .) +  
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10,face='bold'),
          strip.background = element_blank()) +
    labs(title=i,x = "", y = "Gene Expression") +
    theme(plot.title = element_text(hjust = 0.5))+ 
    scale_fill_brewer(palette = "Set1")
  ggsave(p,file=paste0("2.2.Vlnplot_B_",i,".pdf"),width=7,height=4)
}



