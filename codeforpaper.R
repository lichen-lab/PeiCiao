library(dplyr)
library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(cowplot)


##############
# read in data
readSC<-function(samplename){
  data <- Read10X(data.dir = paste(samplename,"/outs/filtered_feature_bc_matrix/",sep=''))
  data <- CreateSeuratObject(counts = data, project = samplename,min.cells = 3,min.features = 100)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
  VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  data<- subset(data, subset = nFeature_RNA > 100  &  nFeature_RNA < 6000 & percent.mt<5)
  data

}


files=list.files('./')
files=files[grep('^D',files)]
ifnb.list=list()
for(ifile in 1:length(files)){
  message(files[ifile])
  data=readSC(files[ifile])
  data@meta.data$condition=files[ifile]
  ifnb.list[files[ifile]]=data
}



#####################################
# perform clustering analysis for D5
ifnb.list.D5=ifnb.list[grep('D5',names(ifnb.list))]


ifnb.list.D5 <- lapply(X = ifnb.list.D5, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


anchors <- FindIntegrationAnchors(ifnb.list.D5, dims = 1:30)
integrated <- IntegrateData(anchors, dims = 1:30)


integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "UMAP")
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated,resolution=0.8)



pdf(file='D5.fig1.pdf',height=4,width=8)
p1 <- DimPlot(integrated, reduction = "UMAP", group.by = "orig.ident")
p2 <- DimPlot(integrated, reduction = "UMAP", label = TRUE)
plot_grid(p1, p2)
dev.off()


pdf(file='D5.fig2.pdf',height=4,width=8)
DimPlot(integrated, reduction = "UMAP", split.by = "orig.ident")
dev.off()



##############################
# perform trajectory analysis

ifnb.list.trace=ifnb.list[grep('3C2I',names(ifnb.list))]

features <- SelectIntegrationFeatures(ifnb.list.trace)

for (i in seq_along(along.with = ifnb.list.trace)) {
  message(i)
  ifnb.list.trace[[i]] <- ScaleData(ifnb.list.trace[[i]], features = features) %>% RunPCA(features = features)
}


anchors <- FindIntegrationAnchors(ifnb.list.trace, reference = c(2), reduction = "rpca", dims = 1:30)
integrated <- IntegrateData(anchors, dims = 1:30)


integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "UMAP")
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated,resolution=0.8)


cds <- as.cell_data_set(integrated)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
unique(colData(cds)$orig.ident)
colData(cds)$time=colData(cds)$orig.ident
colData(cds)$time=factor(colData(cds)$time)
cds <- order_cells(cds)


pdf(file='moncoletrace.days.all.pdf')
plot_cells(cds,
           color_cells_by = "time",
           label_cell_groups=FALSE,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=1.5)
dev.off()



pdf(file='moncoletrace.pseoduotime.all.pdf')
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=1.5)
dev.off()






































