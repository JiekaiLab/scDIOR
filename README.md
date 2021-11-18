

# scDIOR

scDIOR: Single cell RNA-seq  Data IO softwaRe

<br>

<h2 id="0">
    Directory
</h2>

* scDIOR
   * [Overview](#1)
   * [Installing scDIOR](#2)
      * [1. Docker image](#2.1)
      * [2. Conda environment](#2.2)
      * [Version control](#2.3)
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

## Overview  [![top](Figures/top.jpg)](#0)

scDIOR software contains two modules, [dior]() for R and [diopy]() for Python. The data conversion was implemented by a ‘.h5’ file of [HDF5](https://www.hdfgroup.org/) format, which harmonizes the different data types between R and Python. The different aspects of single-cell information were stored in HDF5 group with dataset. scDIOR creates 8 HDF5 groups to store core single-cell information, including data, layers, obs, var, dimR, graphs, uns and spatial.   

![overview](Figures/overview.jpg)

<div id="2"></div>

## Installing scDIOR[![top](Figures/top.jpg)](#0)

Users install and  operate scDIOR following two ways:

1. Docker images are available on the [jiekailab/scdior-image](https://hub.docker.com/r/jiekailab/scdior-image).
2. The environment is created by `conda create` in which scDIOR is installed.



### 1. Docker image

It is recommend to perform scDIOR in docker image, which ensures that the operating environment remains stable. scDIOR image is available on the [jiekailab/scdior-image](https://hub.docker.com/r/jiekailab/scdior-image).

**Brief description**

1. We first built the basic jupyter image which based on [jupyter/base-notebook](https://github.com/jupyter/docker-stacks) (jupyter managing Python and R) and [fixuid](https://github.com/boxboat/fixuid) (fixing user/group mapping issues in containers). This basic image is on [jiekailab/scdior-image:base-jupyter-notebook1.0](https://hub.docker.com/r/jiekailab/scdior-image/tags).
2. Based on our customized basic image, we built scDIOR image again by `Dockerfile`. For the content of `Dockerfile`, it is at this [link](https://github.com/JiekaiLab/scDIOR/blob/main/Dockerfile/Dockerfile).

The current latest image contains the following main analysis platforms and software: 

| R                    | version | Python  | version |
| :------------------- | ------- | ------- | ------- |
| R                    | 4.0.5   | Python  | 3.8.8   |
| Seurat               | 4.0.2   | Scanpy  | 1.8.1   |
| SingleCellExperiment | 1.12.0  | scvelo  | 0.2.3   |
| monocle3             | 1.0.0   | anndata | 0.7.6   |
| dior                 | 0.1.5   | diopy   | 0.5.2   |

### 2. Conda environment

The environment is created by `conda create` in which dior and diopy are installed.

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

 At present, scDIOR is widely compatible with Seurat (v3\~v4) and Scanpy (1.4\~1.8) in different docker image. We configured multiple version docker image (https://hub.docker.com/repository/docker/jiekailab/scdior-image) to confirm that scDIOR can work well between multiple versions of Scanpy and Seurat.

[Demo link]()

| Platform | Software | Version | data IO                 |
| -------- | -------- | ------- | ----------------------- |
| R        | Seurat   | v3~v4   | :ballot_box_with_check: |
| Python   | Scanpy   | 1.4~1.8 | :ballot_box_with_check: |

____



##  Operation Environment 

### 1. Docker image

Pulling the scDIOR image by `docker pull jiekailab/scdior-image`.

1. **Local computer**

   1. Starting the container of scDIOR image

   ```sh
   docker run -p 8840:8888 \ # Port mapping (local computer and docker image)
           --name scdior \
           -e JUPYTER_ENABLE_LAB=yes \
           -e JUPYTER_TOKEN=1234 \
           -v ~/:/home/jovyan \ # Directory mapping
           --rm \
           -it jiekailab/scdior-image
   ```

   2. Opening jupyter notebook in user's  browser 

   ```shell
   localhost:8840
   ```

2. **Remote server **

   1.  Logging in  server through `ssh -L`

   ```shell
   ssh -L localhost:8840:localhost:8853 user@remote_ip
   # localhost:8840:localhost:8853 Port mapping (local and remote server)
   #  user@remote_ip remote useer
   ```

   2. Starting container of scDIOR image

   ```shell
   IMG=jiekailab/scdior-image
   PORT=8853
   PROJECT=scdior
   MEMORY=64g 
   CWD=$(docker inspect $IMG | grep WorkingDir | head -n 1 | sed 's/.* "//;s/"//g;s/,//g') 
   docker run -p $PORT:8888 \ # Port mapping (remote server  and docker image)
           --name $PROJECT \
           -m $MEMORY \
           -u $(id -u):$(id -g) \ # id mapping (user id and user group)
           -e JUPYTER_ENABLE_LAB=yes \
           -e JUPYTER_TOKEN=1234 \
           -v $PWD:$CWD \ # # Directory mapping
           --rm \
           -it $IMG
   ```

   3. Starting jupyter in user's  browser 

   ```shell
   localhost:8840
   ```



### 2. Conda environment

The environment is  activated by `conda activate`.

1. **Local computer**

   1. Activating conda environment

   ```shell
   conda activate conda_env
   ```

   2. Operating  `jupyter lab` or `jupyter notebook` 

   ```shell
   jupyter lab --port=8840
   # or jupyter notebook --port=8840
   ```

   3. Opening jupyter notebook in user's  browser 

   ```shell
   localhost:8840
   ```

2. **Remote server **

   1.  Logging in  server through `ssh -L`

   ```shell
   ssh -L localhost:8840:localhost:8853 user@remote_ip
   # localhost:8840:localhost:8853 Port mapping (local and remote server)
   #  user@remote_ip remote useer
   ```

   2. Activating conda environment

   ```shell
   conda activate conda_env
   ```

   3. Operating  `jupyter lab` or `jupyter notebook` 

   ```shell
   jupyter lab --port=8853
   # or jupyter notebook --port=8853
   ```

   4. Opening jupyter notebook in user's  browser 

   ```shell
   localhost:8840
   ```



## scDIOR demo[![top](Figures/top.jpg)](#0)

Here, we list several demos to show the powerful performance of scDIOR.



### 1. Comparison of trajectory inferences

Users can perform trajectory analysis using Monocle3 in R, then transform the single-cell data to Scanpy in Python using scDIOR, such as expression profiles of spliced and unspliced, as well as cell layout. The expression profile can be used to run dynamical RNA velocity analysis and results can be projected on the layout of Monocle3.

[Demo link](https://github.com/JiekaiLab/scDIOR/tree/main/scdior_demo/Seurat4.0.5_Scanpy1.8.1/1.trajectory_inference)

![1.trajectory_inference](Figures/1.trajectory_inference.png)

### 2. Data IO for batch correction

User can employ single-cell data preprocess and normalization method provided by Scanpy, and utilize batches correction method provided by Seurat.

[Demo link](https://github.com/JiekaiLab/scDIOR/tree/main/scdior_demo/Seurat4.0.5_Scanpy1.8.1/2.batch_correction)

![batch_correct](Figures/batch_correct.png)

### 3. Data IO for spatial omics data

scDIOR supports spatial omics data IO between R and Python platforms.

[Demo link](https://github.com/JiekaiLab/scDIOR/tree/main/scdior_demo/Seurat4.0.5_Scanpy1.8.1/3.spatial_analysis)

![sptail_summary](Figures/sptail_summary.png)

### 4. Extended function

1. The function to load ‘.rds’ file in Python directly;
2. The function to load ‘.h5ad’ file in R directly;
3. Command line  

[Demo link](https://github.com/JiekaiLab/scDIOR/tree/main/scdior_demo/Seurat4.0.5_Scanpy1.8.1/4.scDIOR_extended_function)

![extend_function](Figures/extend_function.png)



## Reference websites [![top](Figures/top.jpg)](#0)

1. jupyter docker stacks: 
   1. https://github.com/jupyter/docker-stacks
   2. https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html
2. fixuid: https://github.com/boxboat/fixuid
3. Seurat: https://satijalab.org/seurat/index.html
4. Scanpy: https://scanpy.readthedocs.io/en/stable/index.html
5. Scvelo: https://scanpy.readthedocs.io/en/stable/index.html
