# scDIOR
scDIOR: Single cell data IO softwaRe

<br>

<h2 id="0">
    Directory
</h2>

* scDIOR
   * [Overview](#1)
   * [Preparation](#2)
      * [Operating Environment](#2.1)
      * [Version control](#2.2)
   * [Getting started](#3)
      * [Starting environment (for docker image)](#3.1)
      * [R loading packages](#3.2)
      * [Python loading packages](#3.3)
   * [Example A](#4)
      * [Loading data with `scvelo` in `Python`](#4.1)
      * [Saving data with `diopy` in `Python`](#4.2)
      * [Loading data with`dior` in `R`](#4.3)
      * [Saving data with`dior` in `R`](#4.4)
      * [Loading `cds_trajectory.h5` with `diopy` in `Python`](#4.5)
   * [Example B](#5)
      * [Loading data with `diopy`in `Python`](#5.1)
      * [Saving data with `diopy` in `Python`](#5.2)
      * [Loading data with`dior` in `R`](#5.3)
   * [Example C](#6)
      * [Loading data from 10X Genomics Spatial Datasets](#6.1)
      * [Saving data by `diopy` in `Python`](#6.2)
      * [Loading data by `dior` in `R`](#6.3)
      * [Saving data by `dior` in `R`](#6.4)
   * [scDIOR extended function](#7)
      * [scDIOR read h5ad file](#7.1)
      * [scDIOR read rds file](#7.2)
      * [scDIOR command line](#7.3)
   * [The scripts link of dior and diopy](#8)
   * [Reference websites](#9)

___

<div id="1"></div>
## Overview [![top](Figures/top.jpg)](#0)

scDIOR software contains two modules, [dior]() for R and [diopy]() for Python. The data transformation was implemented by a ‘.h5’ file of [HDF5](https://www.hdfgroup.org/) format, which harmonizes the different data types between R and Python. The different aspects of single-cell information were stored in HDF5 group with dataset. scDIOR creates 8 HDF5 groups to store core single-cell information, including data, layers, obs, var, dimR, graphs, uns and spatial.   

![overview](Figures/overview.jpg)

<div id="2"></div>
## Preparation [![top](Figures/top.jpg)](#0)

<div id="2.1"></div>
### Operating Environment

**1. Docker image (recommended) :**

It is recommend to perform scDIOR in docker image, which ensures that the operating environment remains stable. scDIOR image is already available for download at Docker Hub, [click here](https://hub.docker.com/repository/docker/jiekailab/scdior-image).  

We Building scDIOR image:

1. we first built the basic jupyter image which based on [jupyter/base-notebook](https://github.com/jupyter/docker-stacks) (jupyter managing Python and R) and [fixuid](https://github.com/boxboat/fixuid) (fixing user/group mapping issues in containers). This basic image is already available for download at Docker Hub, click here.
2.  Based on our customized basic image, we construct scDIOR image again by `Dockerfile`. For the content of `Dockerfile`, please check the file under the dockerfile.

The current latest image contains the following main analysis platforms and software: 

| R                    | version | Python  | version |
| :------------------- | ------- | ------- | ------- |
| R                    | 4.0.5   | Python  | 3.8.8   |
| Seurat               | 4.0.2   | Scanpy  | 1.      |
| SingleCellExperiment | 1.12.0  | scvelo  | 0.2.3   |
| monocle3             | 1.0.0   | anndata | 0.7.6   |
| dior                 | 0.1.4   | diopy   | 0.5.2   |

**2. conda environment (not recommended) :**

The conda environment is changeable without Docker image, so it is not recommended to install the package in this way. The R and Python environments was built by conda, then dior and diopy are installed in R and Python respectively. 

```shell
conda create -n conda_env python=3.8 R=4.0
```

1. R installation:

```R
# in R
install.packages('devtools')
devtools::install_github('JiekaiLab/dior')
# or devtools::install_github('JiekaiLab/dior@HEAD')
```

2. Python installation:

```shell
# in python
pip install diopy
```

<div id="2.2"></div>
### Version control

 At present, scDIOR is widely compatible with Seurat (v3\~v4) and Scanpy (1.4\~1.8) in different docker image. We configured mutitple version docker image (https://hub.docker.com/repository/docker/jiekailab/scdior-image) to confirm that scDIOR can work well between multiple versions of Scanpy and Seurat.dad ag

| Platform | Software | Version | data IO                 |
| -------- | -------- | ------- | ----------------------- |
| R        | Seurat   | v3~v4   | :ballot_box_with_check: |
| Python   | Scanpy   | 1.4~1.8 | :ballot_box_with_check: |

____

<div id="3"></div>
## Getting started [![top](Figures/top.jpg)](#0)

Here, we list the three specific examples and the extended function to show the powerful performance of scDIOR.

* The three examples: 
  1. One can perform trajectory analysis using Monocle3 in R, then transform the single-cell data to Scanpy in Python using scDIOR, such as expression profiles of spliced and unspliced, as well as cell layout. The expression profile can be used to run dynamical RNA velocity analysis and results can be projected on the layout of Monocle3.
  2. One can employ single-cell data preprocess and normalization method provided by Scanpy, and utilize batches correction method provided by Seurat.
  3. scDIOR supports spatial omics data IO between R and Python platforms.

* The extended function:
  1. the function to load ‘.rds’ file in Python directly;
  2. the function to load ‘.h5ad’ file in R directly;
  3. command line  

<div id="3.1"></div>
### Starting environment (for docker image)

1. **Remote server **

   1. Logining server through `ssh L`

   ```shell
   ssh -L localhost:port1:localhost:port2 user@remote_ip
   # port1: local port
   # prot2: remote port
   # user: remote sever user id
   # remote_ip: remote severip
   ```

   2. Starting container of scDIOR image

   ```shell
   IMG=hjfeng/scdior_image:1.2 # 1. 指定使用的镜像
   PORT=port2 # port2
   PROJECT=scdior
   MEMORY=64g 
   CWD=$(docker inspect $IMG | grep WorkingDir | head -n 1 | sed 's/.* "//;s/"//g;s/,//g') 
   docker run -p $PORT:8888 \
           --name $PROJECT \
           -m $MEMORY \
           -u $(id -u):$(id -g) \
           -e JUPYTER_ENABLE_LAB=yes \
           -e JUPYTER_TOKEN=1234 \
           -v $PWD:$CWD \
           --rm \
           -it $IMG
   ```

   3. Starting jupyter in user's  browser 

   ```shell
   localhost:port1
   ```

2. **Local computer**

   1. Staring the container of scDIOR image

   ```shell
   IMG=hjfeng/scdior_image:1.2 # 1. 指定使用的镜像
   PORT=port1# port2
   PROJECT=scdior
   MEMORY=64g 
   CWD=$(docker inspect $IMG | grep WorkingDir | head -n 1 | sed 's/.* "//;s/"//g;s/,//g') 
   docker run -p $PORT:8888 \
           --name $PROJECT \
           -m $MEMORY \
           -u $(id -u):$(id -g) \
           -e JUPYTER_ENABLE_LAB=yes \
           -e JUPYTER_TOKEN=1234 \
           -v $PWD:$CWD \
           --rm \
           -it $IMG
   ```

   2. Starting jupyter in user's  browser 

   ```shell
   localhost:port1
   ```

<div id="3.2"></div>
### R loading packages


```R
# in R
library(Seurat)
library(SingleCellExperiment)
library(dior)
library(monocle3)
library(ggplot2)
```

<div id="3.3"></div>
### Python loading packages

```python
# in python
import scipy
import scanpy as sc
import pandas as pd
import numpy as np
import scvelo as scv
import diopy
```

___

<div id="4"></div>
## Example A [![top](Figures/top.jpg)](#0)

One can perform trajectory analysis using Monocle3 in R, then transform the single-cell data to Scanpy in Python using scDIOR, such as expression profiles of spliced and unspliced, as well as cell layout. The expression profile can be used to run dynamical RNA velocity analysis and results can be projected on the layout of Monocle3.

<div id="4.1"></div>
### Loading data with `scvelo` in `Python`

This data is curated by the `scvelo`, loaded by the code:

```shell
# in python
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

<div id="4.2"></div>
### Saving data with `diopy` in `Python`

```Python
# in python
diopy.output.write_h5(adata = adata, 
                      file = './result/py_write_h5/data_write_velocity.h5')
```

<div id="4.3"></div>
### Loading data with`dior` in `R`

```R
# in R 
sce <- read_h5(file= '.result/py_write_h5/data_write.h5', 
              target.object = 'singlecellexperiment')
cds  <- new_cell_data_set(sce@assays@data@listData$X,
                         cell_metadata = colData(sce),
                         gene_metadata = rowData(sce))
```

* Constructing single-cell trajectories by `monocle3`. More details are available at [monocle3](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/)

1. Pre-process the data

```R
# in R 
cds <- preprocess_cds(cds, 
                      num_dim = 50)
```

2.  Dimensionality reduction 

```R
# in R 
cds <- reduce_dimension(cds)
```

3. Clustering  the cells 

```R
# in R 
cds <- cluster_cells(cds, 
                     cluster_method= 'leiden')
```

4. Learning the trajectory graph and visualization

```R
# in R 
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

5. Adding the Dimension

```R
# in R 
reducedDim(sce,'PCA') <- reducedDim(cds, 'PCA')
reducedDim(sce,'UMAP') <- reducedDim(cds, 'UMAP')
reducedDimNames(sce)<- c('pca','umap','monocle_PCA','monocle_UMAP')
```

<div id="4.4"></div>
### Saving data with`dior` in `R`

```R
# in R 
write_h5(data = sce, 
         file = './result/r_monocle3_result/cds_trajectory.h5', 
         assay.name = 'RNA' ,
         object.type = 'singlecellexperiment')
```

<div id="4.5"></div>
### Loading `cds_trajectory.h5` with `diopy` in `Python`

```python
# in python 
adata = diopy.input.read_h5(file = './result/r_monocle3_result/cds_trajectory.h5')
```

* RNA Velocity analysis. More details are available at [scvelo](https://scvelo.readthedocs.io/VelocityBasics/)

1. Preprocess the Data

```python
# in python 
scv.pp.filter_and_normalize(adata, 
                            min_shared_counts=20, 
                            n_top_genes=2000)
scv.pp.moments(adata, 
               n_pcs=30, 
               n_neighbors=30)
```

2. Estimate RNA velocity

```python
# in python 
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
```

3. Project the velocities

```Python
# in python 
scv.pl.velocity_embedding_stream(adata, 
                                 basis='monocle_UMAP')
```

![trajectory_inference_by_scvelo](Figures/trajectory_inference_by_scvelo.png)

___

<div id="5"></div>
## Example B [![top](Figures/top.jpg)](#0)

One can employ single-cell data preprocess and normalization method provided by Scanpy, and utilize batches correction method provided by Seurat.

<div id="5.1"></div>
### Loading data with `diopy`in `Python`

The data is curated by `scanpy`, loaded by the code:

```python
# in python 
adata = sc.read('data/pancreas.h5ad',
                backup_url='https://www.dropbox.com/s/qj1jlm9w10wmt0u/pancreas.h5ad?dl=1')
```

1. Visualization

```python
# in python 
sc.pl.umap(adata, 
           color=['batch', 'celltype'], 
           palette=sc.pl.palettes.vega_20_scanpy)
```

![batch_data_in_scanpy](Figures/batch_data_in_scanpy.png)

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

![batch_data_corrected_by_seurat](Figures/batch_data_corrected_by_seurat.jpg)

___

<div id="6"></div>
## Example C [![top](Figures/top.jpg)](#0)

<div id="6.1"></div>
### Loading data from 10X Genomics Spatial Datasets

Downloading the spatial data set from [10X Genomics Spatial Datasets](https://support.10xgenomics.com/spatial-gene-expression/datasets)

* Analysis and visualization of spatial transcriptomics data by [scanpy spatial](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html)

1. Reading the data

```python
# in python
adata = sc.read_visium('./data/V1_Human_Lymph_Node')
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
```

2. QC and preprocessing

```python
# in python
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20]
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)
```

3.  Normalize Visium counts data 

```python
# in python
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
```

4. Manifold embedding and clustering 

```python
# in python
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")
```

5. Visualization in spatial coordinates

```python
# in python
sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2"], save='.spatial_py_cluster_gene.png')
```

![spatail_analysis_by_scanpy](Figures/spatail_analysis_by_scanpy.png)

<div id="6.2"></div>
### Saving data with`diopy` in `Python`

```python
# in python
diopy.output.write_h5(adata=adata, 
                      file='./result/spatial_data_scanpy.h5',
                      assay_name='spatial')
```

<div id="6.3"></div>
### Loading data with `dior` in `R`

```R
# in R
adata = read_h5(file = './result/spatial_data_scanpy.h5', 
                assay.name = 'spatial')
```

1. Visualization in spatial coordinates

```R
# in R
options(repr.plot.width=8, repr.plot.height=8)
SpatialDimPlot(adata, 
               label = TRUE, 
               label.size = 3, 
               group.by='clusters',  
               pt.size.factor = 1)
```

![r_spatial_cluster](Figures/r_spatial_cluster.jpg)

```R
# in R
SpatialFeaturePlot(adata, 
                   features = c('CR2'), 
                   pt.size.factor = 1)
```

![r_spatial_gene_cr2](Figures/r_spatial_gene_cr2.jpg)

<div id="6.4"></div>
### Saving data with `dior` in `R`

```R
# in R
write_h5(adata, 
         file = './result/spatial_data_scanpy_v2.h5', 
         assay.name = 'spatial')
```

load the data by diopy in Python

```python
# in python
adata = diopy.input.read_h5(file = './result/spatial_data_scanpy_v2.h5',
                            assay_name='spatial')
```

___

<div id="7"></div>
## scDIOR extended function [![top](Figures/top.jpg)](#0)

<div id="7.1"></div>
### scDIOR read h5ad file

Reading the h5ad file in R. `dior::read_h5ad` function will create a file with `_tmp.h5` suffix.

```R
# ~/.conda/envs/vev1/bin/R
adata = read_h5ad(file = './data/data_test_batch.h5ad', 
                  target.object = 'seurat', 
                  assay_name = 'RNA')
```

<div id="7.2"></div>
### scDIOR read rds file

Reading the rds file in Python. `diopy.input.read_rds` function will create a file with `_tmp.h5` suffix.

```python
# in python
adata = diopy.input.read_rds(file = './result/r_monocle3_result/sce_trajectory.rds',
                             object_type='singlecellexperiment',
                             assay_name='RNA')
adata
# AnnData object with n_obs × n_vars = 3696 × 27998
#     obs: 'clusters_coarse', 'clusters', 'S_score', 'G2M_score'
#     var: 'highly_variable_genes'
#     uns: 'clusters_colors'
#     obsm: 'X_mono_pca', 'X_mono_umap', 'X_pca', 'X_umap'
#     layers: 'spliced', 'unspliced'
```

<div id="7.3"></div>
### scDIOR command line

 ScDIOR uses the command line to convert different data by calling `scdior`.  

`usage: scdior [-h] -i INPUT -o OUTPUT -t TARGET -a ASSAY_NAME`

`-i,--input` The existing filename for different platforms, such as rds (R) or h5ad (Python)

`-o,--output`  The filename that needs to be converted, such as from rds to h5ad or from h5ad to rds

`-t,--target` The target object for R, such as seruat or singlecellexperiment

`-a,--assay_name` The primary data types, such as scRNA data or spatial data

* Example scdior convert the h5ad to the rds

```shell
# in shell
$ scdior -i ./data_test_batch.h5ad -o ./data_test_batch.rds -t seurat -a RNA
# ...loading the h5ad file...
# Warning: No columnames present in cell embeddings, setting to 'PCA_1:50'
# Warning: No columnames present in cell embeddings, setting to 'UMAP_1:2'
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

# ...saving the rds file...
# ...complete....
```

```R
# in R
library(Seurat)
adata <- readRDS('./data_test_batch.rds')
adata
# An object of class Seurat
# 24516 features across 14693 samples within 1 assay
# Active assay: RNA (24516 features, 0 variable features)
#  2 dimensional reductions calculated: pca, umap
```

* Example scdior convert the rds to the h5ad

```shell
# in shell
$ scdior -i ./data_test_batch.rds -o ./data_test_batch.h5ad -t seurat -a RNA
# ...loading the rds file...
# ...saving the h5ad file...
# ...complete....
```

```python
# in python
import scanpy as sc
adata = sc.read('./data_test_batch.h5ad')
adata
# AnnData object with n_obs × n_vars = 14693 × 24516
#     obs: 'celltype', 'sample', 'n_genes', 'batch', 'n_counts', 'louvain'
#     var: 'n_cells-0', 'n_cells-1', 'n_cells-2', 'n_cells-3'
#     obsm: 'X_pca', 'X_umap'
#     obsp: 'connectivities', 'distances'
```

___

<div id="8"></div>
## The scripts link of dior and diopy [![top](Figures/top.jpg)](#0)

1. 


# djaihga g


<div id="9"></div>
## Reference websites [![top](Figures/top.jpg)](#0)

1. jupyter docker stacks: 
   1. https://github.com/jupyter/docker-stacks
   2. https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html
2. fixuid: https://github.com/boxboat/fixuid
3. Seurat: https://satijalab.org/seurat/index.html
4. Scanpy: https://scanpy.readthedocs.io/en/stable/index.html
5. Scvelo: https://scanpy.readthedocs.io/en/stable/index.html
