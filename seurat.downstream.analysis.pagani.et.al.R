library(Seurat)

outs_folder <- 'path/to/outputs'
sets <- list()
set.day = c('00', '07', '42', '42')
for (i in c('Sample_352-AH-1', 'Sample_352-AH-2', 'Sample_352-AH-3', 'Sample_352-AH-4')){
    raw_data <- Read10X(data.dir =
                          paste0('path/to/sample/', i, '-3GEX/filtered_feature_bc_matrix'))
    colnames(x = raw_data) <- paste(paste(i, sep = '.'), colnames(x = raw_data), sep = '_') 
    print(dim(raw_data))
    sets[i] <- CreateSeuratObject(counts = raw_data, project = "Pagani.et.al", min.cells = 1, min.features = 0)
    sets[i][["perc.mt"]] <- PercentageFeatureSet(sets[[sample.tags[i]]], pattern = "^mt-")
    sets[i][["logUMI"]] <- log(sets[[sample.tags[i]]][["nCount_RNA"]])
    sets[i][[day]] = set.day[i]
}

rm(raw_data)

max_genes = c(5000, 7500, 6000, 6000)
for (i in 1:length(sets)){
  sets[[i]] <- subset(sets[[i]],
                      subset = nFeature_RNA < max_genes[i] & nFeature_RNA > 500 & sets[[i]][["perc.mt"]] < 25)
  print(dim(sets[[i]]))
  }

first.letter.up <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

s.genes <- first.letter.up(cc.genes$s.genes)
g2m.genes <- first.letter.up(cc.genes$g2m.genes)

for (i in 1:length(sets)) {
  sets[[i]] <- NormalizeData(sets[[i]], verbose = FALSE)
  sets[[i]] <- CellCycleScoring(sets[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  sets[[i]] <- FindVariableFeatures(sets[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
  }

sets.anchors <- FindIntegrationAnchors(object.list = sets, dims = 1:50)
sets <- IntegrateData(anchorset = sets.anchors, dims = 1:50)

DefaultAssay(sets) <- 'integrated'

sets <- ScaleData(sets, vars.to.regress = c('S.Score', 'G2M.Score', 'perc.mt'), verbose = TRUE)

sets <- RunPCA(sets, verbose = FALSE, npcs = 150)
ElbowPlot(sets, ndims = 150)

sets <- RunUMAP(sets, dims = 1:75, umap.method = 'umap-learn', metric = 'correlation')

sets <- FindNeighbors(sets, dims = 1:75)
sets <- FindClusters(sets, res = 0.35, n.start = 100, n.iter = 100, algorithm = 1)

sets$TdTomatopos <- as.factor(sets[['RNA']]@counts['TdTomato',] > 0)
