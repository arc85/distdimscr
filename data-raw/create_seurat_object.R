## Create Seurat object for dimdistscr
## Dec 13 2021

## Load Seurat
library(Seurat)

## Read in samples and create metdata
files.to.read <- grep("\\.R",list.files("data-raw"),invert=TRUE,value=TRUE)
ser.list <- vector("list",length=length(files.to.read))

for (i in 1:length(files.to.read)) {

  tmp.sample <- Read10X(paste("data-raw",files.to.read[i],sep="/"))
  tmp.meta <- data.frame(sample_id=rep(NA,ncol(tmp.sample)),sample_type=rep(NA,ncol(tmp.sample)))
  tmp.meta$sample_id <- rep(files.to.read[i],length(tmp.meta$sample_id))
  tmp.meta$sample_type <- ifelse(grepl("PBMC",tmp.meta$sample_id),"PBMC","Tonsil")
  rownames(tmp.meta) <- colnames(tmp.sample)

  ser.list[[i]] <- CreateSeuratObject(tmp.sample,meta.data=tmp.meta)

}

## Merge into one Seurat Object
ser <- merge(ser.list[[1]],ser.list[2:length(ser.list)])

## Seruat clustering and visualization workflow
ser <- NormalizeData(ser)
ser <- FindVariableFeatures(ser)
ser <- ScaleData(ser)
ser <- RunPCA(ser)
ElbowPlot(ser)

ser <- RunUMAP(ser,dims=1:10)
ser <- FindNeighbors(ser,dims=1:10)
ser <- FindClusters(ser,res=0.3)

## Identify cell types
DimPlot(ser,group.by="sample_type")
DimPlot(ser,label=T)
FeaturePlot(ser,c("CD3D","CD8A","CD4","CD14","FCGR3A","MS4A1","MZB1","IL3RA","CLEC10A","HBB"))

## Add cell type metadata
data.to.add <- vector("logical",length=ncol(ser))

data.to.add[ser@meta.data$RNA_snn_res.0.3=="0" |
              ser@meta.data$RNA_snn_res.0.3=="7" |
              ser@meta.data$RNA_snn_res.0.3=="6"] <- "B cells"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="8"] <- "Plasmablasts"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="1" |
              ser@meta.data$RNA_snn_res.0.3=="3" |
              ser@meta.data$RNA_snn_res.0.3=="2" |
              ser@meta.data$RNA_snn_res.0.3=="5" ] <- "NK and T cells"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="4"] <- "CD14 monocytes"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="10"] <- "pDCs"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="11"] <- "RBCs"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="9"] <- "CD1C DCs"

ser[["cell_types"]] <- data.to.add

## Visualize cell types
DimPlot(ser,group.by="cell_types",split.by="sample_type")

## Subcluster to identify NK and T cells
ser.t <- ser[,ser@meta.data$cell_types=="NK and T cells"]
ser.t <- FindVariableFeatures(ser.t)
ser.t <- ScaleData(ser.t)
ser.t <- RunPCA(ser.t)
ElbowPlot(ser.t)

ser.t <- RunUMAP(ser.t,dims=1:10)
ser.t <- FindNeighbors(ser.t,dims=1:10)
ser.t <- FindClusters(ser.t,res=0.5)

DimPlot(ser.t,label=T)
FeaturePlot(ser.t,c("CD3D","CD8A","CD4","FCGR3A"))

## Add T cell and NK cell types
data.to.add <- vector("logical",length=ncol(ser.t))

data.to.add[ser.t@meta.data$RNA_snn_res.0.5=="1" |
              ser.t@meta.data$RNA_snn_res.0.5=="2" |
              ser.t@meta.data$RNA_snn_res.0.5=="11" |
              ser.t@meta.data$RNA_snn_res.0.5=="7" |
              ser.t@meta.data$RNA_snn_res.0.5=="5" |
              ser.t@meta.data$RNA_snn_res.0.5=="0" |
              ser.t@meta.data$RNA_snn_res.0.5=="3" |
              ser.t@meta.data$RNA_snn_res.0.5=="6"] <- "CD4 T cells"

data.to.add[ser.t@meta.data$RNA_snn_res.0.5=="4" |
              ser.t@meta.data$RNA_snn_res.0.5=="9" |
              ser.t@meta.data$RNA_snn_res.0.5=="12" |
              ser.t@meta.data$RNA_snn_res.0.5=="10"] <- "CD8 T cells"

data.to.add[ser.t@meta.data$RNA_snn_res.0.5=="8"] <- "NK cells"

ser.t[["cell_types"]] <- data.to.add

## Check out T cell types
DimPlot(ser.t,group.by="cell_types")

## Add back to overall object
ser@meta.data$cell_types[match(colnames(ser.t),colnames(ser))] <- ser.t@meta.data$cell_types
DimPlot(ser,group.by="cell_types",split.by="sample_type")

## Metadata
overall.metadata <- ser@meta.data
format(object.size(overall.metadata),unit="MB")

## UMAP
overall.umap <- Embeddings(ser,reduction="umap")
format(object.size(overall.umap),unit="MB")

## PCA embeddings
overall.pca <- Embeddings(ser,reduction="pca")[,1:10]
format(object.size(overall.pca),unit="MB")

## Save minimum files
save(overall.metadata,file="data/overall_metadata.RData")
save(overall.umap,file="data/overall_umap.RData")
save(overall.pca,file="data/overall_pca.RData")
