# Data IO for batch correction

Users can employ single-cell data preprocess and normalization method provided by Scanpy, and utilize batches correction method provided by Seurat.



## Links



https://github.com/fenghuijian/fenghuijian.github.io/tree/master/doc/scdior_demo/Seurat4.0.5_Scanpy1.8.1/2.batch_correction





## Demo

### Loading data with `diopy`in `Python`

The data is curated by `scanpy`, loaded by the code:

```python
# in python 
adata = sc.read('data/pancreas.h5ad',backup_url='https://www.dropbox.com/s/qj1jlm9w10wmt0u/pancreas.h5ad?dl=1')
```

1. Visualization

```python
# in python 
sc.pl.umap(adata, 
           color=['batch', 'celltype'], 
           palette=sc.pl.palettes.vega_20_scanpy)
```

![batch_data_in_scanpy](../../../Figures/batch_data_in_scanpy.png)

<div id="5.2"></div>

### Saving data with `diopy` in `Python`

```python
# in python 
diopy.output.write_h5(adata, 
                      file = './result/batch_effect_data.h5', 
                      save_X=False) # Select not to save adata_all.X, because that's scale data,
```

<div id="5.3"></div>

### Loading data with `dior` in `R`

```R
# in R
adata <- read_h5(file = './result/batch_effect_data.h5',
                 assay.name = 'RNA', 
                 target.object = 'seurat')
adata@meta.data$batch <- as.character(adata@meta.data$batch)
```

* Batch effect corrected by Seurat protocol. More details are available at [Seurat](https://satijalab.org/seurat/articles/get_started.html)

1. Dataset preprocessing: Splitting the combined object into a list

```R
# in R
adata_list <- SplitObject(adata, 
                          split.by = "batch")
```

2. Dataset preprocessing:  Variable feature selection based on a variance stabilizing transformation (`"vst"`) 

```R
# in R
adata_list <- lapply(X = adata_list, FUN = function(x) {
    # x <- NormalizeData(x) The data is normal data and does not need to be normalized
    x <- FindVariableFeatures(x, 
                              selection.method = "vst", 
                              nfeatures = 2000)
})
```

3. Select integration features

```R
# in R
features <- SelectIntegrationFeatures(object.list = adata_list)
```

4. Integration of cell datasets

```R
# in R
adata_anchors <- FindIntegrationAnchors(object.list = adata_list,
                                     anchor.features = features)
adata_combined <- IntegrateData(anchorset = adata_anchors)
```

5.  Downstream analysis for Integration data

```R
# in R
DefaultAssay(adata_combined) <- "integrated"
adata_combined <- ScaleData(adata_combined, verbose = FALSE)
adata_combined <- RunPCA(adata_combined, npcs = 30, verbose = FALSE)
adata_combined <- RunUMAP(adata_combined, reduction = "pca", dims = 1:30)
```

6.  Visualization 

```R
# in R 
options(repr.plot.width=20, repr.plot.height=8)
DimPlot(adata_combined, reduction = "umap", group.by = c("batch", 'celltype'))
```

![batch_data_corrected_by_seurat](../../../Figures/batch_data_corrected_by_seurat.jpg)