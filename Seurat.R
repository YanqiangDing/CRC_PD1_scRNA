#################################################
# Single cell RNA sequencing
# Seurat
# 2023
#################################################

library("Rcpp")
library(dplyr)
library(Seurat)
library(cowplot)
library(patchwork)
require(ggplot2)
library(dendextend)
require(ComplexHeatmap)
require(circlize)

method="tsne"


dir.create("./2_results")
dir.create("./2_results/1_Data_preprocessing")
dir.create("./2_results/2_Linear_dimensionality_reduction_analysis")
dir.create("./2_results/3_Non_linear_dimensionality_reduction_analysis")
dir.create("./2_results/4_cluster_biomarkers")
dir.create("./2_results/5_cluster_DEgenes")
dir.create("./2_results/6_oncogenic_genes/")


###------------------ALL samples------------------------------------------------------------------------

PA.data <- Read10X(data.dir = "/storage/yanqiangding/20220210_scRNA_BoYi_Jianning/2_count_cellranger/from_com/MC38_PA/filtered_feature_bc_matrix")
dim(PA.data)

CON.data <- Read10X(data.dir = "/storage/yanqiangding/20220210_scRNA_BoYi_Jianning/2_count_cellranger/from_com/MC38_CON/filtered_feature_bc_matrix")
dim(CON.data)

## Create Seurat Object
PA <- CreateSeuratObject(counts = PA.data, project = "PA", min.cells = 3, min.features = 200)

CON <- CreateSeuratObject(counts = CON.data, project = "CON", min.cells = 3, min.features = 200)


CON=readRDS(file = "../filter-lm_1018/CON-filter-D.rds")   #"CON"
PA=readRDS(file = "../filter-lm_1018/PA-filter-D.rds")     #"PA"


##
CON[["percent.mt"]] <- PercentageFeatureSet(CON, pattern = "^mt-")   # genes started with "MT-" is mitochondria genes.
p=VlnPlot(CON, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0)
ggsave("./2_results/1_Data_preprocessing/CON_VlnPlot.jpg", plot = p, width = 8, height = 5)

CON <- subset(CON, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
p=VlnPlot(CON, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0)
ggsave("./2_results/1_Data_preprocessing/CON_VlnPlot_subset.jpg", plot = p, width = 8, height = 5)


##
PA[["percent.mt"]] <- PercentageFeatureSet(PA, pattern = "^mt-")   # genes started with "MT-" is mitochondria genes.
p=VlnPlot(PA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0)
ggsave("./2_results/1_Data_preprocessing/PA_VlnPlot.jpg", plot = p, width = 8, height = 5)

PA <- subset(PA, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
p=VlnPlot(PA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0)
ggsave("./2_results/1_Data_preprocessing/PA_VlnPlot_subset.jpg", plot = p, width = 8, height = 5)


sclist <- c(CON, PA)
for (i in 1:length(sclist)) {

	sclist[[i]] <- NormalizeData(sclist[[i]], verbose = FALSE)
	sclist[[i]] <- FindVariableFeatures(sclist[[i]], selection.method = "vst", 
					    nfeatures = 2000, verbose = FALSE)
}

#--k.filter = 20
SelectIntegrationFeatures(object.list = sclist)
#FindIntegrationAnchors
sc.anchors <- FindIntegrationAnchors(object.list = sclist, dims = 1:20, k.anchor = 5, k.filter = 20)
#IntegrateData
sc.integrated <- IntegrateData(anchorset = sc.anchors, dims = 1:20)

sc.integrated <- ScaleData(sc.integrated, verbose = TRUE)

sc.integrated <- RunPCA(sc.integrated, verbose = TRUE)
# Plot PCA
p=DimPlot(sc.integrated, reduction = "pca",split.by = 'ident')
ggsave("./2_results/2_Linear_dimensionality_reduction_analysis/PCA.jpg", plot = p, width = 8, height = 5)


sc.integrated <- JackStraw(sc.integrated, num.replicate = 100)   # need about 15min
sc.integrated <- ScoreJackStraw(sc.integrated, dims = 1:20)

pdf(file = "./2_results/2_Linear_dimensionality_reduction_analysis/JackStraw_dim20_k20.pdf", page="a4")
JackStrawPlot(sc.integrated, dims = 1:20)
dev.off()

pdf(file = "./2_results/2_Linear_dimensionality_reduction_analysis/Elbow_dim20_k20.pdf", page="a4")
ElbowPlot(sc.integrated)
dev.off()


sc.integrated <- FindNeighbors(sc.integrated, dims = 1:20)
sc.integrated <- FindClusters(sc.integrated, resolution = 0.5)

# dimensionality_reduction
sc.integrated <- RunUMAP(sc.integrated, dims = 1:20)
sc.integrated <- RunTSNE(sc.integrated, dims = 1:20)
p=DimPlot(sc.integrated, reduction = "umap",label = TRUE, split.by="orig.ident")
ggsave("./2_results/4_cluster_biomarkers/uamp_split.jpg", plot = p, width = 8, height = 5)

p=DimPlot(sc.integrated, reduction = "tsne",label = TRUE, split.by="orig.ident")
ggsave("./2_results/4_cluster_biomarkers/tsne_split.jpg", plot = p, width = 8, height = 5)


# find markers for every cluster compared to all remaining cells, report only the positive ones
sc.integrated.markers <- FindAllMarkers(sc.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(file = "./2_results/4_cluster_biomarkers/Markers.csv", sc.integrated.markers)

sc.integrated.markers_top5 = sc.integrated.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(file = "./2_results/4_cluster_biomarkers/Markers_top5.csv", sc.integrated.markers_top5)

#sc.integrated.markers_top20 = sc.integrated.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(file = "./2_results/4_cluster_biomarkers/Markers_top20.csv", sc.integrated.markers_top20)

### marker genes heatmap
genes = c("Gsn", "Dcn", "C3", "Pi16", "Col3a1", "Serping1", "lgfbp6", "Mfap5", "Lum", "Ccl2", "Tnfaip6", "Sod3", "Ccl21a", "Fabb4", "Ctla2a", "Plvap", "Lyve1", "Egfl7", "Emcn", "Pecam1", "Lrg1", "Aqp1", "Mmrn1", "Nts", "Ecscr", "Esam", "Gng11", "Krt14", "Krt5", "Lgals7", "Krt17", "S100a8", "Krt15", "Fabp5", "Krt6a", "Gsta4", "Krt13", "Mt4", "Perp", "Gpx2", "Cd3e", "cd3d", "cd5", "cd28", "cd4", "cd8a", "Trbc1", "Gzma", "Klrb1c", "Ncr1", "Klra4", "Eomes", "Gzmb", "Fasl", "Klrk1", "Cd19", "Cd79a", "Cd79b", "Ms4a1", "Fcmr", "Pax5", "Cd22", "Ly6g", "S100a8", "S100a9", "Mmp9", "Csf3r", "Mmp8", "Cxcr2", "Ly6c2", "Fcgr1", "Csf1r", "Zeb2", "Sirpa", "Ccr2", "Xcr1", "Clec9a", "Itgax", "Irf8", "clec10a", "Ly6c2", "Irf8", "clec10a", "Zeb2", "Ccr9", "siglech", "Bst2", "Tcf4", "Pacsin1", "Mrc1", "Mertk", "cd8a", "Trbc1", "Gzma", "Gzmb", "cd4", "Trbc1", "Trbc1", "Gzma", "Klrb1c", "Ncr1", "Gzmb", "Gzma", "Klrb1c", "Ncr1", "Klra4", "Eomes", "Gzmb", "Fasl", "Klrk1", "Foxp3", "cd4", "Trbc1", "Arg2", "Il1β")

pdf(file = "./2_results/4_cluster_biomarkers/DoHeatmap_marker.pdf", page="a4")
DoHeatmap(sc.integrated, features = genes) + NoLegend()
dev.off()

pdf(file = "./2_results/4_cluster_biomarkers/DoHeatmap_marker_group.pdf", page="a4")
DoHeatmap(sc.integrated, features = genes, group.by = "ident",) + NoLegend()
dev.off()

exp=subset(sc.integrated, features = genes)

# metadata
head(exp@meta.data)

# expression data
exp[["RNA"]][1:3,]

#Calculating the average gene expression within a cluster
cluster.averages <- AverageExpression(sc.integrated)
head(cluster.averages[["RNA"]][, 1:5])

pdf(file = "./2_results/4_cluster_biomarkers/DoHeatmap_marker_exp.pdf", page="a4")
DoHeatmap(cluster.averages, features = genes, size = 3,
              draw.lines = FALSE)
dev.off()



gene_exp=cluster.averages[["RNA"]][row.names(cluster.averages[["RNA"]]) %in% genes,]
write.csv(file = "./2_results/4_cluster_biomarkers/marker_gene_exp.csv", gene_exp)

cd45=cluster.averages[["RNA"]][row.names(cluster.averages[["RNA"]]) %in% "Ptprc",]
write.csv(file = "./2_results/4_cluster_biomarkers/marker_cd45_Ptprc_gene_exp.csv", cd45)

write.csv(file = "./2_results/4_cluster_biomarkers/all_marker_gene_exp.csv", cluster.averages[["RNA"]])


## cluster cell number
cluster_cell_number = table(sc.integrated$seurat_clusters)
write.csv(file = "./2_results/4_cluster_biomarkers/cluster_cell_number.csv", cluster_cell_number)


# save sc.integrated$seurat_clusters
write.csv(file = "./2_results/4_cluster_biomarkers/seurat_clusters.csv", sc.integrated$seurat_clusters)

write.csv(file = "./2_results/4_cluster_biomarkers/orig.ident.csv", sc.integrated$orig.ident)

write.csv(file = "./2_results/4_cluster_biomarkers/nCount_RNA.csv", sc.integrated$nCount_RNA)

write.csv(file = "./2_results/4_cluster_biomarkers/nFeature_RNA.csv", sc.integrated$nFeature_RNA)

write.csv(file = "./2_results/4_cluster_biomarkers/percent.mt.csv", sc.integrated$percent.mt)

write.csv(file = "./2_results/4_cluster_biomarkers/integrated_snn_res.0.5.csv", sc.integrated$integrated_snn_res.0.5)


sc.integrated_CON=subset(sc.integrated, orig.ident == "CON")

sc.integrated_PA=subset(sc.integrated, orig.ident == "PA")


CON_Cluster=table(sc.integrated_CON$seurat_clusters)

PA_Cluster=table(sc.integrated_PA$seurat_clusters)


write.csv(file = "./2_results/4_cluster_biomarkers/CON_Cluster.csv", CON_Cluster)

write.csv(file = "./2_results/4_cluster_biomarkers/PA_Cluster.csv", PA_Cluster)


## 
saveRDS(sc.integrated, file = './sc.integrated_DimRed.Rds')



#sc.integrated = readRDS(file = "./sc.integrated_DimRed.Rds")



##---------------- Non immune cell--------------------------------

summary(GetAssayData(sc.integrated, assay = "RNA")["Ptprc",])

GetAssayData(sc.integrated, assay = "RNA")["Ptprc", GetAssayData(sc.integrated, assay = "RNA")["Ptprc",] == 0]

count_RNA=GetAssayData(sc.integrated, assay = "RNA")

class(count_RNA)

count_RNA[ , count_RNA["Ptprc", ] == 0]
dim(count_RNA[ , count_RNA["Ptprc", ] == 0])

NOimm_cell_id=colnames(count_RNA[ , count_RNA["Ptprc", ] == 0])

sc.integrated_NOimm <- sc.integrated[ ,NOimm_cell_id]

saveRDS(sc.integrated_NOimm, file = './sc.integrated_NOimm_DimRed.Rds')


##-----------immune cell------------------------------------------------------------------------

dir.create("./2_results/immune_cell/")
dir.create("./2_results/immune_cell/1_Data_preprocessing")
dir.create("./2_results/immune_cell/2_Linear_dimensionality_reduction_analysis")
dir.create("./2_results/immune_cell/3_Non_linear_dimensionality_reduction_analysis")
dir.create("./2_results/immune_cell/4_cluster_biomarkers")
dir.create("./2_results/immune_cell/5_cluster_DEgenes")

dim_num=20

imm_cell_id=colnames(count_RNA[ , count_RNA["Ptprc", ] > 0])

sc.integrated_imm <- sc.integrated[ ,imm_cell_id]

saveRDS(sc.integrated_imm, file = './sc.integrated_imm_DimRed.Rds')


sc.integrated_imm <- RunPCA(sc.integrated_imm, verbose = TRUE)
# Plot PCA
p=DimPlot(sc.integrated_imm, reduction = "pca",split.by = 'ident')
ggsave("./2_results/immune_cell/2_Linear_dimensionality_reduction_analysis/PCA.jpg", plot = p, width = 8, height = 5)


sc.integrated_imm <- JackStraw(sc.integrated_imm, num.replicate = 100)  
sc.integrated_imm <- ScoreJackStraw(sc.integrated_imm, dims = 1:dim_num)

pdf(file = "./2_results/immune_cell/2_Linear_dimensionality_reduction_analysis/JackStraw.pdf", page="a4")
JackStrawPlot(sc.integrated_imm, dims = 1:dim_num)
dev.off()

pdf(file = "./2_results/immune_cell/2_Linear_dimensionality_reduction_analysis/Elbow.pdf", page="a4")
ElbowPlot(sc.integrated_imm)
dev.off()


sc.integrated_imm <- FindNeighbors(sc.integrated_imm, dims = 1:dim_num)
sc.integrated_imm <- FindClusters(sc.integrated_imm, resolution = 0.5)


# dimensionality_reduction
sc.integrated_imm <- RunUMAP(sc.integrated_imm, dims = 1:dim_num)
sc.integrated_imm <- RunTSNE(sc.integrated_imm, dims = 1:dim_num)


p=DimPlot(sc.integrated_imm, reduction = "umap",label = TRUE, split.by="orig.ident")
ggsave("./2_results/immune_cell/4_cluster_biomarkers/uamp_split.jpg", plot = p, width = 8, height = 5)

p=DimPlot(sc.integrated_imm, reduction = "tsne",label = TRUE, split.by="orig.ident")
ggsave("./2_results/immune_cell/4_cluster_biomarkers/tsne_split.jpg", plot = p, width = 8, height = 5)


# find markers for every cluster compared to all remaining cells, report only the positive ones
sc.integrated_imm.markers <- FindAllMarkers(sc.integrated_imm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/Markers.csv", sc.integrated_imm.markers)

sc.integrated_imm.markers_top5 = sc.integrated_imm.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/Markers_top5.csv", sc.integrated_imm.markers_top5)

saveRDS(sc.integrated_imm, file = './2_results/immune_cell/sc.integrated_imm_DimRed.Rds')

p=DimPlot(sc.integrated_imm, reduction = "tsne",
	  label = FALSE, label.size = 2,
	  split.by = "orig.ident",
	  pt.size = 0.05) +
	NoLegend()
ggsave("./2_results/immune_cell/4_cluster_biomarkers/tsne_NOlabel_split_celltype3.jpg", plot = p, width = 8, height = 5)

pdf(file = "./2_results/immune_cell/4_cluster_biomarkers/tsne_NOlabel_split_celltype3.pdf", width = 8, height = 5)
DimPlot(sc.integrated_imm, reduction = "tsne",
          label = FALSE, label.size = 2,
          split.by = "orig.ident",
          pt.size = 0.05) +
        NoLegend()

dev.off()


p=DimPlot(sc.integrated_imm, reduction = "tsne",
           label = TRUE, label.size = 2,
           split.by = "orig.ident",
           pt.size = 0.05) +
        NoLegend()
ggsave("./2_results/immune_cell/4_cluster_biomarkers/tsne_label_split_celltype2.jpg", plot = p, width = 8, height = 5)

pdf(file = "./2_results/immune_cell/4_cluster_biomarkers/tsne_label_split_celltype2.pdf", width = 8, height = 5)
DimPlot(sc.integrated_imm, reduction = "tsne",
	 label = TRUE, label.size = 2,
	 split.by = "orig.ident",
	 pt.size = 0.05) +
        NoLegend()

dev.off()


p=DimPlot(sc.integrated_imm, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("./2_results/immune_cell/4_cluster_biomarkers/tsne_label_NOsplit_celltype.jpg", plot = p, width = 8, height = 5)

pdf(file = "./2_results/immune_cell/4_cluster_biomarkers/tsne_label_NOsplit_celltype.pdf", page="a4")
DimPlot(sc.integrated_imm, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()


p=DimPlot(sc.integrated_imm, reduction = "tsne", label = TRUE, split.by = "orig.ident", pt.size = 0.5) + NoLegend()
ggsave("./2_results/immune_cell/4_cluster_biomarkers/tsne_label_split_celltype.jpg", plot = p, width = 8, height = 5)

pdf(file = "./2_results/immune_cell/4_cluster_biomarkers/tsne_label_split_celltype.pdf", page="a4")
DimPlot(sc.integrated_imm, reduction = "tsne", label = TRUE, split.by = "orig.ident", pt.size = 0.5) + NoLegend()
dev.off()

p=DimPlot(sc.integrated_imm, reduction = "umap", label = TRUE, split.by = "orig.ident", pt.size = 0.5) + NoLegend()
ggsave("./2_results/immune_cell/4_cluster_biomarkers/umap_label_split_celltype.jpg", plot = p, width = 8, height = 5)

pdf(file = "./2_results/immune_cell/4_cluster_biomarkers/umap_label_split_celltype.pdf", page="a4")
DimPlot(sc.integrated_imm, reduction = "umap", label = TRUE, split.by = "orig.ident", pt.size = 0.5) + NoLegend()
dev.off()


### marker genes EXPRESSION

genes = c("Gsn", "Dcn", "C3", "Pi16", "Col3a1", "Serping1", "lgfbp6", "Mfap5", "Lum", "Ccl2", "Tnfaip6", "Sod3", 
          "Ccl21a", "Fabb4", "Ctla2a", "Plvap", "Lyve1", "Egfl7", "Emcn", "Pecam1", "Lrg1", "Aqp1", "Mmrn1", "Nts", 
          "Ecscr", "Esam", "Gng11", "Krt14", "Krt5", "Lgals7", "Krt17", "S100a8", "Krt15", "Fabp5", "Krt6a", "Gsta4", 
          "Krt13", "Mt4", "Perp", "Gpx2", "Cd3e", "Cd3d", "Cd5", "Cd28", "Cd4", "Cd8a", "Trbc1", "Gzma", "Klrb1c", "Ncr1", 
          "Klra4", "Eomes", "Gzmb", "Fasl", "Klrk1", "Cd19", "Cd79a", "Cd79b", "Ms4a1", "Fcmr", "Pax5", "Cd22", "Ly6g", 
          "S100a8", "S100a9", "Mmp9", "Csf3r", "Mmp8", "Cxcr2", "Ly6c2", "Fcgr1", "Csf1r", "Zeb2", "Sirpa", "Ccr2", 
          "Xcr1", "Clec9a", "Itgax", "Irf8", "clec10a", "Ly6c2", "Irf8", "clec10a", "Zeb2", "Ccr9", "siglech", "Bst2", 
          "Tcf4", "Pacsin1", "Mrc1", "Mertk", "Ncr1", "Klrb1c", "Ncr1", "Klra4", "Eomes", "Fasl", "Klrk1", "Foxp3", 
          "Arg2", "Il1β", "Il1b", "Ifng", "Krt86", "Tnfrsf10b", "Vegfa", "Il1b", "Cd209a", "Camp", "Itgam", "Adgre1", "Siglech")

pdf(file = "./2_results/immune_cell/4_cluster_biomarkers/DoHeatmap_marker.pdf", page="a4")
DoHeatmap(sc.integrated_imm, features = genes) + NoLegend()
dev.off()

pdf(file = "./2_results/immune_cell/4_cluster_biomarkers/DoHeatmap_marker_group.pdf", page="a4")
DoHeatmap(sc.integrated_imm, features = genes, group.by = "ident",) + NoLegend()
dev.off()

exp=subset(sc.integrated_imm, features = genes)

# metadata
head(exp@meta.data)

# expression data
exp[["RNA"]][1:3,]

#Calculating the average gene expression within a cluster
cluster.averages <- AverageExpression(sc.integrated_imm)
head(cluster.averages[["RNA"]][, 1:5])


gene_exp=cluster.averages[["RNA"]][row.names(cluster.averages[["RNA"]]) %in% genes,]
write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/marker_gene_exp.csv", gene_exp)

gene_exp=cluster.averages[["RNA"]][row.names(cluster.averages[["RNA"]]), ]
write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/marker_ALL_gene_exp.csv", gene_exp)


## cluster cell number
cluster_cell_number = table(sc.integrated_imm$seurat_clusters)
write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/cluster_cell_number.csv", cluster_cell_number)


# save sc.integrated_imm$seurat_clusters
write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/seurat_clusters.csv", sc.integrated_imm$seurat_clusters)

write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/orig.ident.csv", sc.integrated_imm$orig.ident)

write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/nCount_RNA.csv", sc.integrated_imm$nCount_RNA)

write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/nFeature_RNA.csv", sc.integrated_imm$nFeature_RNA)

write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/percent.mt.csv", sc.integrated_imm$percent.mt)

write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/integrated_snn_res.0.5.csv", sc.integrated_imm$integrated_snn_res.0.5)


sc.integrated_imm_CON=subset(sc.integrated_imm, orig.ident == "CON")

sc.integrated_imm_PA=subset(sc.integrated_imm, orig.ident == "PA")


CON_Cluster=table(sc.integrated_imm_CON$seurat_clusters)

PA_Cluster=table(sc.integrated_imm_PA$seurat_clusters)


write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/CON_Cluster.csv", CON_Cluster)

write.csv(file = "./2_results/immune_cell/4_cluster_biomarkers/PA_Cluster.csv", PA_Cluster)
##---------



#-----------T cell----------------------------------------------------------------------------------------------------------------------------------------------------

dir.create("./2_results/Tcell/")
dir.create("./2_results/Tcell/1_Data_preprocessing")
dir.create("./2_results/Tcell/2_Linear_dimensionality_reduction_analysis")
dir.create("./2_results/Tcell/3_Non_linear_dimensionality_reduction_analysis")
dir.create("./2_results/Tcell/4_cluster_biomarkers")
dir.create("./2_results/Tcell/5_cluster_DEgenes")

dim_num=20

t_cell_id=colnames(count_RNA[ , count_RNA["Trbc1", ] > 0])


sc.integrated_t <- sc.integrated[ ,t_cell_id]

saveRDS(sc.integrated_t, file = './sc.integrated_t_DimRed.Rds')

sc.integrated_t <- readRDS(file = './sc.integrated_t_DimRed.Rds')

sc.integrated_t <- RunPCA(sc.integrated_t, verbose = TRUE)
# Plot PCA
p=DimPlot(sc.integrated_t, reduction = "pca",split.by = 'ident')
ggsave("./2_results/Tcell/2_Linear_dimensionality_reduction_analysis/PCA.jpg", plot = p, width = 8, height = 5)


sc.integrated_t <- JackStraw(sc.integrated_t, num.replicate = 100)   # need about 15min
sc.integrated_t <- ScoreJackStraw(sc.integrated_t, dims = 1:dim_num)

pdf(file = "./2_results/Tcell/2_Linear_dimensionality_reduction_analysis/JackStraw.pdf", page="a4")
JackStrawPlot(sc.integrated_t, dims = 1:dim_num)
dev.off()

pdf(file = "./2_results/Tcell/2_Linear_dimensionality_reduction_analysis/Elbow.pdf", page="a4")
ElbowPlot(sc.integrated_t)
dev.off()


sc.integrated_t <- FindNeighbors(sc.integrated_t, dims = 1:dim_num)
sc.integrated_t <- FindClusters(sc.integrated_t, resolution = 0.5)


# dimensionality_reduction
sc.integrated_t <- RunUMAP(sc.integrated_t, dims = 1:dim_num)
sc.integrated_t <- RunTSNE(sc.integrated_t, dims = 1:dim_num)


p=DimPlot(sc.integrated_t, reduction = "umap",label = TRUE, split.by="orig.ident")
ggsave("./2_results/Tcell/4_cluster_biomarkers/uamp_split.jpg", plot = p, width = 8, height = 5)

p=DimPlot(sc.integrated_t, reduction = "tsne",label = TRUE, split.by="orig.ident")
ggsave("./2_results/Tcell/4_cluster_biomarkers/tsne_split.jpg", plot = p, width = 8, height = 5)


# find markers for every cluster compared to all remaining cells, report only the positive ones
sc.integrated_t.markers <- FindAllMarkers(sc.integrated_t, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/Markers.csv", sc.integrated_t.markers)

sc.integrated_t.markers_top5 = sc.integrated_t.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/Markers_top5.csv", sc.integrated_t.markers_top5)


saveRDS(sc.integrated_t, file = './2_results/Tcell/sc.integrated_t_DimRed.Rds')



p=DimPlot(sc.integrated_t, reduction = "tsne",
	  label = FALSE, label.size = 2,
	  split.by = "orig.ident",
	  pt.size = 0.05) +
	NoLegend()
ggsave("./2_results/Tcell/4_cluster_biomarkers/tsne_NOlabel_split_celltype3.jpg", plot = p, width = 8, height = 5)

pdf(file = "./2_results/Tcell/4_cluster_biomarkers/tsne_NOlabel_split_celltype3.pdf", width = 8, height = 5)
DimPlot(sc.integrated_t, reduction = "tsne",
          label = FALSE, label.size = 2,
          split.by = "orig.ident",
          pt.size = 0.05) +
        NoLegend()

dev.off()


p=DimPlot(sc.integrated_t, reduction = "tsne",
           label = TRUE, label.size = 2,
           split.by = "orig.ident",
           pt.size = 0.05) +
        NoLegend()
ggsave("./2_results/Tcell/4_cluster_biomarkers/tsne_label_split_celltype2.jpg", plot = p, width = 8, height = 5)

pdf(file = "./2_results/Tcell/4_cluster_biomarkers/tsne_label_split_celltype2.pdf", width = 8, height = 5)
DimPlot(sc.integrated_t, reduction = "tsne",
	 label = TRUE, label.size = 2,
	 split.by = "orig.ident",
	 pt.size = 0.05) +
        NoLegend()

dev.off()


p=DimPlot(sc.integrated_t, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("./2_results/Tcell/4_cluster_biomarkers/tsne_label_NOsplit_celltype.jpg", plot = p, width = 8, height = 5)

pdf(file = "./2_results/Tcell/4_cluster_biomarkers/tsne_label_NOsplit_celltype.pdf", page="a4")
DimPlot(sc.integrated_t, reduction = "tsne", label = TRUE,  pt.size = 0.5) + NoLegend()
dev.off()


p=DimPlot(sc.integrated_t, reduction = "tsne", label = TRUE, split.by = "orig.ident", pt.size = 0.5) + NoLegend()
ggsave("./2_results/Tcell/4_cluster_biomarkers/tsne_label_split_celltype.jpg", plot = p, width = 8, height = 5)

pdf(file = "./2_results/Tcell/4_cluster_biomarkers/tsne_label_split_celltype.pdf", page="a4")
DimPlot(sc.integrated_t, reduction = "tsne", label = TRUE, split.by = "orig.ident", pt.size = 0.5) + NoLegend()
dev.off()

p=DimPlot(sc.integrated_t, reduction = "umap", label = TRUE, split.by = "orig.ident", pt.size = 0.5) + NoLegend()
ggsave("./2_results/Tcell/4_cluster_biomarkers/umap_label_split_celltype.jpg", plot = p, width = 8, height = 5)

pdf(file = "./2_results/Tcell/4_cluster_biomarkers/umap_label_split_celltype.pdf", page="a4")
DimPlot(sc.integrated_t, reduction = "umap", label = TRUE, split.by = "orig.ident", pt.size = 0.5) + NoLegend()
dev.off()


### marker genes EXPRESSION

genes1 = c("Cd3e", "Cd3d", "Cd5", "Cd28", "Cd4", "Cd8a", "Trbc1", "Mki67", "Stmn1", "Pclaf", "Pcna", "Xcl1", "Ncr1", "Cd7", "Cxr6", "Ncr1", "Klrb1c", "Ly6g", "S100a8", "S100a9", "Cd19", "Cd79a", "Cd79b", "Ms4a1", "Xcr1", "Clec9a", "Itgax", "Irf8", "Itgax", "Clec10a", "Ly6c2", "Irf8", "Clec10a", "Zeb2", "Csf1r", "Spn", "Sirpa", "Zeb2", "Ly6c2", "Fcgr1", "Csf1r", "Zeb2", "Sirpa", "Mafb", "Fcgr1", "Adgre1", "Clec4f", "Timd4", "Zeb2", "Sirpa", "Csf1r", "Cd5l", "Pdcd1", "Cd274", "Mrc1", "Mertk", "Foxp3", "Gzma", "Gzmb", "Cd8b1", "Ptprc")

genes = c("Gsn", "Dcn", "C3", "Pi16", "Col3a1", "Serping1", "lgfbp6", "Mfap5", "Lum", "Ccl2", "Tnfaip6", "Sod3", "Ccl21a", "Fabb4", "Ctla2a", "Plvap", "Lyve1", "Egfl7", "Emcn", "Pecam1", "Lrg1", "Aqp1", "Mmrn1", "Nts", "Ecscr", "Esam", "Gng11", "Krt14", "Krt5", "Lgals7", "Krt17", "S100a8", "Krt15", "Fabp5", "Krt6a", "Gsta4", "Krt13", "Mt4", "Perp", "Gpx2", "Cd3e", "cd3d", "cd5", "cd28", "cd4", "cd8a", "Trbc1", "Gzma", "Klrb1c", "Ncr1", "Klra4", "Eomes", "Gzmb", "Fasl", "Klrk1", "Cd19", "Cd79a", "Cd79b", "Ms4a1", "Fcmr", "Pax5", "Cd22", "Ly6g", "S100a8", "S100a9", "Mmp9", "Csf3r", "Mmp8", "Cxcr2", "Ly6c2", "Fcgr1", "Csf1r", "Zeb2", "Sirpa", "Ccr2", "Xcr1", "Clec9a", "Itgax", "Irf8", "clec10a", "Ly6c2", "Irf8", "clec10a", "Zeb2", "Ccr9", "siglech", "Bst2", "Tcf4", "Pacsin1", "Mrc1", "Mertk", "cd8a", "Trbc1", "Gzma", "Gzmb", "cd4", "Trbc1", "Trbc1", "Gzma", "Klrb1c", "Ncr1", "Gzmb", "Gzma", "Klrb1c", "Ncr1", "Klra4", "Eomes", "Gzmb", "Fasl", "Klrk1", "Foxp3", "cd4", "Trbc1", "Arg2", "Il1β")

pdf(file = "./2_results/Tcell/4_cluster_biomarkers/DoHeatmap_marker.pdf", page="a4")
DoHeatmap(sc.integrated_t, features = genes) + NoLegend()
dev.off()

pdf(file = "./2_results/Tcell/4_cluster_biomarkers/DoHeatmap_marker_group.pdf", page="a4")
DoHeatmap(sc.integrated_t, features = genes, group.by = "ident",) + NoLegend()
dev.off()

exp=subset(sc.integrated_t, features = genes)

# metadata
head(exp@meta.data)

# expression data
exp[["RNA"]][1:3,]

#Calculating the average gene expression within a cluster
cluster.averages <- AverageExpression(sc.integrated_t)
head(cluster.averages[["RNA"]][, 1:5])


gene_exp=cluster.averages[["RNA"]][row.names(cluster.averages[["RNA"]]) %in% genes,]
write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/marker_gene_exp.csv", gene_exp)


gene_exp=cluster.averages[["RNA"]][row.names(cluster.averages[["RNA"]]),]
write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/all_marker_gene_exp.csv", gene_exp)


## cluster cell number
cluster_cell_number = table(sc.integrated_t$seurat_clusters)
write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/cluster_cell_number.csv", cluster_cell_number)

# save sc.integrated_t$seurat_clusters
write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/seurat_clusters.csv", sc.integrated_t$seurat_clusters)

write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/orig.ident.csv", sc.integrated_t$orig.ident)

write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/nCount_RNA.csv", sc.integrated_t$nCount_RNA)

write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/nFeature_RNA.csv", sc.integrated_t$nFeature_RNA)

write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/percent.mt.csv", sc.integrated_t$percent.mt)

write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/integrated_snn_res.0.5.csv", sc.integrated_t$integrated_snn_res.0.5)


sc.integrated_t_CON=subset(sc.integrated_t, orig.ident == "CON")

sc.integrated_t_PA=subset(sc.integrated_t, orig.ident == "PA")

CON_Cluster=table(sc.integrated_t_CON$seurat_clusters)

PA_Cluster=table(sc.integrated_t_PA$seurat_clusters)


write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/CON_Cluster.csv", CON_Cluster)

write.csv(file = "./2_results/Tcell/4_cluster_biomarkers/PA_Cluster.csv", PA_Cluster)