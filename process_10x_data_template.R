# on rstudio profiler on gdp
# on rstudio gdp
# cd /scratch/Oncology_Shared/oncology_scRNA
# scp -r nbatada@rbiodev.gilead.com:/resbioinfo/data/Oncology/Team_Shared/public_scRNAseq/Zhang2024_CRC_GSE236581 .
DIR='/scratch/Oncology_Shared/oncology_scRNA/Zhang2024_CRC_GSE236581/filtered_feature_bc_matrix'


library(Seurat)
sobj=CreateSeuratObject(counts=Read10X(data.dir=DIR),  min.cells=10, min.features=500) # names.delim = "_" to extract out patient name from barcodes which will be put in orig.ident
#
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
#VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
sobj <- subset(sobj, subset = nFeature_RNA < 800 & percent.mt < 5)
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000) # In Seurat v5, Normalized values are stored in sobj[["RNA"]]$data.
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(sobj), 10)

# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(sobj)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2

## linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1 
#By default, only variable features are scaled.
#You can specify the features argument to scale additional features
# The results of this are stored in sobj[["RNA"]]$scale.data 
all.genes <- rownames(sobj)
sobj <- ScaleData(sobj, features = all.genes)
sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))
sobj <- FindNeighbors(sobj, dims = 1:30) # as in paper
sobj <- FindClusters(sobj, resolution=0.5)

# too slow
#sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE)
#sobj.markers %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1)

meta_file_publication='/scratch/Oncology_Shared/oncology_scRNA/Zhang2024_CRC_GSE236581/GSE236581_CRC-ICB_metadata.txt'
meta.data <- read.csv(meta_file_publication, sep=' ')
meta.data=data.frame(meta.data)
sobj <- AddMetaData(sobj, meta.data)
head(sobj@meta.data,1)

# export all
saveRDS(sobj, file = "/scratch/Oncology_Shared/oncology_scRNA/Zhang2024_CRC_GSE236581/Zhang2024_CRC_GSE236581_full.rds")

# export subset
unique(sobj@meta.data$MajorCellType)
Idents(sobj) <- "MajorCellType" # define current ident from a column of interest in sobj@meta.data
sobj_tcell=subset(sobj, ident=c("T"))
class(sobj_tcell)
saveRDS(sobj_tcell, file = "/scratch/Oncology_Shared/oncology_scRNA/Zhang2024_CRC_GSE236581/Zhang2024_CRC_GSE236581_tcells.rds")

unique(sobj@meta.data$SubCellType)
Idents(sobj) <- "SubCellType" # define current ident from a column of interest in sobj@meta.data
sobj_tregs=subset(sobj, idents=c("c11_CD4_Treg_FOXP3","c12_CD4_Treg_KLRB1","c13_CD4_Treg_TNFRSF9"))
saveRDS(sobj_tregs, file = "/scratch/Oncology_Shared/oncology_scRNA/Zhang2024_CRC_GSE236581/Zhang2024_CRC_GSE236581_tregs.rds")


library(ggplot2)
library(ggpubr) # stat_compare_mean()
#install.packages('ggpubr')
head(sobj_tcell@meta.data,1)
Idents(sobj_tcell) <- "SubCellType" # define current ident from a column of interest in sobj@meta.data

sobj_tcell@meta.data$SubCellType=as.factor(sobj_tcell@meta.data$SubCellType)

plot=VlnPlot(sobj_tcell, features=c("FOXP3"), sort=T, pt.size=.1, layer="data", log=TRUE) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) # +  ggpubr::stat_compare_means() + ylim(-2,15)
ggsave('vlnplot_foxp3.pdf', height=7, width=12, plot=plot)


DotPlot(sobj_tcell, features=c("FOXP3", "CD4", "CD3E", "CCR8"))

VlnPlot(sobj_tcell, features=c("CCR8"), sort=T, pt.size=.1,layer="data", log=TRUE) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) # +  ggpubr::stat_compare_means() + ylim(-2,15)
ggsave('vlnplot_ccr8.pdf', height=7, width=12, plot=plot)

VlnPlot(sobj_tcell, features=c("CD3E"), sort=T, pt.size=.1,layer="data", log=TRUE) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) # +  ggpubr::stat_compare_means() + ylim(-2,15)
VlnPlot(sobj_tcell, features=c("CD4"), sort=T, pt.size=.1,layer="data", log=TRUE) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) # +  ggpubr::stat_compare_means() + ylim(-2,15)
VlnPlot(sobj_tcell, features=c("CD8A"), sort=T, pt.size=.1,layer="data", log=TRUE) + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) # +  ggpubr::stat_compare_means() + ylim(-2,15)


