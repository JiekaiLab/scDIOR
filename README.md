# scDIOR
scDIOR: Single cell data IO softwaRe

[toc]

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

```python
# env/ python
import scipy
import scanpy as sc
import pandas as pd
import numpy as np
import scvelo as scv
import diopy
```



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



### Load the data in R

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



### load the `cds_trajectory.h5` by diopy inPython

```Python
mono = diopy.input.read_h5(file = './result/r_monocle3_result/cds_trajectory.h5')
```

1. filter and normalization and 

```Python
scv.pp.filter_and_normalize(mono, 
                            min_shared_counts=20, 
                            n_top_genes=2000)
scv.pp.moments(mono, 
               n_pcs=30, 
               n_neighbors=30)

scv.tl.velocity(mono)
scv.tl.velocity_graph(mono)
```



```Python
scv.pl.velocity_embedding_stream(mono, 
                                 basis='mono_umap',
                                 save='.scvelo_trajectory.png')
```

![trajectory_inference_by_scvelo](Figures/trajectory_inference_by_scvelo.png)







[dior](https://jiekailab.github.io/scDior/sc_data_IO_r.html)

[diopy](https://jiekailab.github.io/scDior/sc_data_IO_python.html)



