# scDior
scDior: Single cell data IO softwaRe

## installation

### R

```R
install.packages('devtools')
devtools::install_github('JiekaiLab/RIOH5@HEAD')
```

### Python

```shell
git clone https://github.com/JiekaiLab/PyIOH5.git
cd PyIOH5
python setup.py install
```



## Getting started

### 1. The data is saved in R and can be read in R & Python

#### R: library

```R
library(Seurat)
library(SingleCellExperiment)
library(hdf5r)
library(RIOH5)
```

#### R: save data into the h5 file

```R
# read the data from rds file
data = readRDS('./seurat.rds' )
# save the data into the h5 file 
write_h5(data = data, object.type = 'seurat',file = './seurat_r.h5')
```

#### R: read data from the h5 file

```R
# read the h5 file for seurat object
data_seurat = read_h5(target.object = 'seurat', file = './seurat_r.h5')
# read the h5 file for singlcellexperiment object
data_sce = read_h5(target.object = 'singlecellexperiment', file = './seurat_r.h5')
```

 #### Python: read data from the h5 file

```python
# import 
import PyIOH5 as myh5
import scanpy as sc
# read the h5 file for scanpy object
data = myh5.input.read_h5(file='./seurat_r.h5')
```



### 2. The data is saved in Python and can be read in R & Python.

#### Python: import

```Python
import PyIOH5 as myh5
import scanpy as sc
```

#### Python: save data into the h5 file

```Python
# read the data from h5ad file
adata = sc.read( file = './scanpy.h5ad')
# save the data into the h5 file 
myh5.ouput.write_h5(adata=adata, file = './py_scanpy.h5')
```

#### Python: read data from the h5 file

```Python
adata_r = myh5.input.read_h5(file ='./py_scanpy.h5')
```

#### R: read data from the h5 file 

```R
# read the h5 file for seurat object
data_seurat = read_h5(target.object = 'seurat', file = './py_scanpy.h5')
# read the h5 file for singlcellexperiment object
data_sce = read_h5(target.object = 'singlecellexperiment', file = './py_scanpy.h5')
```



 ## See our tutorial page for more details 

**Download the data and tutorial** 

```shell
git clone https://github.com/JiekaiLab/JiekaiLab.github.io.git
cd JiekaiLab.github.io/scDior
```

**Or see the online tutorial**

[RIOH5](https://jiekailab.github.io/scDior/sc_data_IO_r.html)

[PyIOH5](https://jiekailab.github.io/scDior/sc_data_IO_python.html)



