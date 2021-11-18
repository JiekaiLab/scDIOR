## scDIOR extended function



## Link

[4.1.extended_function_in_Python](https://fenghuijian.github.io/doc/scdior_demo/Seurat4.0.5_Scanpy1.8.1/4.scDIOR_extended_function/4.1.extended_function_in_Python.html)

shows that the function of  Pratial extraction  and Loading rds file in Python.  [4.2.extended_function_in_R](https://fenghuijian.github.io/doc/scdior_demo/Seurat4.0.5_Scanpy1.8.1/4.scDIOR_extended_function/4.2.extended_function_in_R.html)

shows that the function of  Pratial extraction and Loading h5ad file in R.



## Demo

### scDIOR read h5ad file

Reading the h5ad file in R. `dior::read_h5ad` function will create a file with `_tmp.h5` suffix.

```R
# in R
adata_seurat = read_h5ad(file = './adata_Python.h5ad', 
                  target.object = 'seurat', 
                  assay_name = 'RNA')
adata_seurat
# An object of class Seurat 
# 83994 features across 100 samples within 3 assays 
# Active assay: RNA (27998 features, 0 variable features)
#  2 other assays present: spliced, unspliced
#  2 dimensional reductions calculated: pca, umap
```

```R
# in R
adata_sce = read_h5ad(file ='./adata_Python.h5ad', 
                         target.object = 'singlecellexperiment',
                         assay_name = 'RNA')
adata_sce
# class: SingleCellExperiment 
# dim: 27998 100 
# metadata(0):
# assays(3): X spliced unspliced
# rownames(27998): Xkr4 Gm37381 ... Gm20837 Erdr1
# rowData names(1): highly_variable_genes
# colnames(100): TTCTACAGTACTCTCC GGACGTCCAGGGAGAG ... CACACCTCATTGTGCA
#   ACTGCTCAGTTACCCA
# colData names(4): clusters_coarse clusters S_score G2M_score
# reducedDimNames(2): pca umap
# altExpNames(0):
```



### scDIOR read rds file

Reading the rds file in Python. `diopy.input.read_rds` function will create a file with `_tmp.h5` suffix.

```python
# in python
adata = diopy.input.read_rds(file = './adata_R.rds',
                             object_type='seurat',
                             assay_name='RNA')
adata
# AnnData object with n_obs × n_vars = 100 × 27998
#     obs: 'clusters_coarse', 'clusters', 'S_score', 'G2M_score'
#     var: 'highly_variable_genes'
#     obsm: 'X_pca', 'X_umap'
#     layers: 'spliced', 'unspliced'
#     obsp: 'distances', 'connectivities'
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
$ scdior -i ./adata_test.h5ad -o ./adata_test.rds -t seurat -a RNA
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
adata <- readRDS('./adata_test.rds')
adata
# An object of class Seurat
# 24516 features across 14693 samples within 1 assay
# Active assay: RNA (24516 features, 0 variable features)
#  2 dimensional reductions calculated: pca, umap
```

* Example scdior convert the rds to the h5ad

```shell
# in shell
$ scdior -i ./adata_test.rds -o ./adata_test.h5ad -t seurat -a RNA
# ...loading the rds file...
# ...saving the h5ad file...
# ...complete....
```

```python
# in python
import scanpy as sc
adata = sc.read('./adata_test.h5ad')
adata
# AnnData object with n_obs × n_vars = 14693 × 24516
#     obs: 'celltype', 'sample', 'n_genes', 'batch', 'n_counts', 'louvain'
#     var: 'n_cells-0', 'n_cells-1', 'n_cells-2', 'n_cells-3'
#     obsm: 'X_pca', 'X_umap'
#     obsp: 'connectivities', 'distances'
```

