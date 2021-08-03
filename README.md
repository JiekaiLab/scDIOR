# scDIOR
scDIOR: Single cell data IO softwaRe

* [scDIOR](#scdior)
   * [overview](#overview)
   * [installation](#installation)
      * [R](#r)
      * [Python](#python)
   * [Getting started](#getting-started)
   * [scDIOR Fig. 3A example](#scdior-fig-3a-example)
      * [load the data in Python](#load-the-data-in-python)
      * [save the data with diopy](#save-the-data-with-diopy)
      * [Load the data by dior in R](#load-the-data-by-dior-in-r)
      * [save the data by dior](#save-the-data-by-dior)
      * [load the cds_trajectory.h5 by diopy in Python](#load-the-cds_trajectoryh5-by-diopy-in-python)
   * [scDIOR Fig. 3B example](#scdior-fig-3b-example)
      * [load the data in Python](#load-the-data-in-python-1)
      * [save the data by diopy in Python](#save-the-data-by-diopy-in-python)
      * [load the data by dior in R](#load-the-data-by-dior-in-r-1)
   * [scDIOR Fig. 3C example](#scdior-fig-3c-example)
      * [load the data from <a href="https://support.10xgenomics.com/spatial-gene-expression/datasets" rel="nofollow">10X Genomics</a>](#load-the-data-from-10x-genomics)
      * [save the data by diopy in Python](#save-the-data-by-diopy-in-python-1)



## overview

![overview](Figures/overview.jpg)

## installation

### R

```R
install.packages('devtools')
devtools::install_github('JiekaiLab/dior')
# devtools::install_github('JiekaiLab/dior@HEAD')
```

### Python

`pip` is recommended

```shell
pip install diopy
```



## Getting started

Python 

```python
# ~/.conda/envs/vev1/bin/python
import scipy
import scanpy as sc
import pandas as pd
import numpy as np
import scvelo as scv
import diopy
```

R

```R
# ~/.conda/envs/vev1/bin/R
library(Seurat)
library(SingleCellExperiment)
library(dior)
library(monocle3)
library(ggplot2)
```



___

____

____



## scDIOR Fig. 3A example

### load the data in Python

This data is curated by the `scvelo`, loaded by the code:

```shell
adata = scv.datasets.pancreas()
adata
# AnnData object with n_obs × n_vars = 3696 × 27998
#     obs: 'clusters_coarse', 'clusters', 'S_score', 'G2M_score'
#     var: 'highly_variable_genes'
#     uns: 'clusters_coarse_colors', 'clusters_colors', 'day_colors', 'neighbors', 'pca'
#     obsm: 'X_pca', 'X_umap'
#     layers: 'spliced', 'unspliced'
#     obsp: 'distances', 'connectivities'
```

### save the data with diopy

```Python
diopy.output.write_h5(adata = adata, 
                      file = '/data1/home/jkchen/hjfeng/Projects/h5/paper_script/result/py_write_h5/data_write_velocity.h5')
```

### Load the data by dior in R

```R
sce = read_h5(file= '/data1/home/jkchen/hjfeng/Projects/h5/paper_script/result/py_write_h5/data_write.h5', 
              target.object = 'singlecellexperiment')
cds  <- new_cell_data_set(sce@assays@data@listData$X,
                         cell_metadata = colData(sce),
                         gene_metadata = rowData(sce))
```

1. Preprocess a cds to prepare for trajectory inference by monocle3 

```R
cds <- preprocess_cds(cds, 
                      num_dim = 50)
```

2. Compute a projection of a cell_data_set object into a lower dimensional space with non-linear dimension reduction methods by monocle3

```R
cds <- reduce_dimension(cds)
```

3. Plotting 

```R
plot_cells(cds, 
           label_groups_by_cluster=FALSE,  
           color_cells_by = "clusters", 
           rasterize = T,
           cell_stroke=0, 
           cell_size = 2,  
           reduction_method = 'UMAP', 
           group_label_size =8)
```

4. Clustering

```R
cds <- cluster_cells(cds, 
                     cluster_method= 'leiden')
```

5. Trajectory inference

```R
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "clusters",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           reduction_method = 'UMAP',
           cell_size = 1, 
           group_label_size =8)
```

![trajectory_inference_by_monocle3](Figures/trajectory_inference_by_monocle3.png)

### save the data by dior

```R
write_h5(data = sce, 
         file = './result/r_monocle3_result/cds_trajectory.h5', assay.name = 'RNA' ,
         object.type = 'singlecellexperiment')
```



### load the `cds_trajectory.h5` by diopy in Python

```python
mono = diopy.input.read_h5(file = './result/r_monocle3_result/cds_trajectory.h5')
```

1. filter and normalization and 

```python
scv.pp.filter_and_normalize(mono, 
                            min_shared_counts=20, 
                            n_top_genes=2000)
scv.pp.moments(mono, 
               n_pcs=30, 
               n_neighbors=30)

scv.tl.velocity(mono)
scv.tl.velocity_graph(mono)
```

2. fi

```Python
scv.pl.velocity_embedding_stream(mono, 
                                 basis='mono_umap',
                                 save='.scvelo_trajectory.png')
```

![trajectory_inference_by_scvelo](Figures/trajectory_inference_by_scvelo.png)

____

____

____



## scDIOR Fig. 3B example

### load the data in Python

```python
adata_all = sc.read('data/pancreas.h5ad', backup_url='https://www.dropbox.com/s/qj1jlm9w10wmt0u/pancreas.h5ad?dl=1')
```

1. ploting

```python
sc.pl.umap(adata_all, 
           color=['batch', 'celltype'], 
           palette=sc.pl.palettes.vega_20_scanpy)
```

![batch_data_in_scanpy](Figures/batch_data_in_scanpy.png)

### save the data by diopy in Python

```python
diopy.output.write_h5(adata_all, 
                      file = './result/batch_effect_data.h5', 
                      save_X=False) # no save adata_all.X because it is the scale data
```

### load the data by dior in R

```R
data_batch <- read_h5(file = './result/batch_effect_data.h5',
                      assay.name = 'RNA', 
                      target.object = 'seurat')

data_batch@meta.data$batch <- as.character(data_batch@meta.data$batch)
```

* Batch effect corrected by Seurat protocol. More details are available at [Seurat](https://satijalab.org/seurat/articles/get_started.html)

1. Dataset preprocessing: Splitting the combined object into a list

```R
db_list <- SplitObject(data_batch, 
                       split.by = "batch")
```

2. Dataset preprocessing:  Variable feature selection based on a variance stabilizing transformation (`"vst"`) 

```R
db_list <- lapply(X = db_list, FUN = function(x) {
    # x <- NormalizeData(x) 该数据已经是normal数据，不需要进行normal
    x <- FindVariableFeatures(x, 
                              selection.method = "vst", 
                              nfeatures = 2000)
})
```

3. 

```R
features <- SelectIntegrationFeatures(object.list = db_list)
```



```R
db_anchors <- FindIntegrationAnchors(object.list = db_list,
                                     anchor.features = features)
```



```R
db_combined <- IntegrateData(anchorset = db_anchors)
```



```R
DefaultAssay(db_combined) <- "integrated"
```



```R
db_combined <- ScaleData(db_combined, verbose = FALSE)
db_combined <- RunPCA(db_combined, npcs = 30, verbose = FALSE)
db_combined <- RunUMAP(db_combined, reduction = "pca", dims = 1:30)
```



```R
options(repr.plot.width=20, repr.plot.height=8)
DimPlot(db_combined, reduction = "umap", group.by = c("batch", 'celltype'))
```

![batch_data_corrected_by_seurat](Figures/batch_data_corrected_by_seurat.jpg)



## scDIOR Fig. 3C example

### load the data from [10X Genomics](https://support.10xgenomics.com/spatial-gene-expression/datasets)

* Analysis and visualization of spatial transcriptomics data by [scanpy spatail](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html)

```python
adata = sc.read_visium('./data/V1_Human_Lymph_Node')
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
```



```pyhon
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20]
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)
```



```python
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
```



```python
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")
```



```python
sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2"], save='.spatial_py_cluster_gene.png')
```

![spatail_analysis_by_scanpy](Figures/spatail_analysis_by_scanpy.png)







### save the data by diopy in Python











[dior](https://jiekailab.github.io/scDior/sc_data_IO_r.html)

[diopy](https://jiekailab.github.io/scDior/sc_data_IO_python.html)



