FROM jiekailab/scdior-image:base-jupyter-notebook1.0
# this base iamge is from docker hub jiekailab/scdior-image:base-jupyter-notebook1.0
# https://hub.docker.com/repository/docker/jiekailab/scdior-image
ARG DEBIAN_FRONTEND=noninteractive 
RUN conda config --add channels bioconda \
    && mamba install -y -c conda-forge \
        git \
        leidenalg \
        numba=0.51.0 \
        pytables \
        python-igraph \
        seaborn \
        scikit-learn \
        statsmodels \
        jupyter \
        jupyterlab \
        r-base=4.0.5 \
        r-remotes \
        r-devtools \
        r-doparallel \
        r-dt \
        r-ggpubr \
        r-hdf5r \
        r-hmisc \
        r-irkernel \
        r-matrix \
        r-pheatmap \
        r-seurat \
        r-getopt \
        r-ggrastr \
        r-loom \
        bioconductor-delayedmatrixstats \
        bioconductor-batchelor \
        bioconductor-limma \
        bioconductor-singlecellexperiment \
        bioconductor-biocversion \
        bioconductor-biocgenerics \
        bioconductor-scran \
        bioconductor-scater \
        bioconductor-loomexperiment \
        r-leidenbase \
        r-monocle3 \
    && pip install \
        scanpy==1.8.1 \
        scvelo \
        diopy \
        loompy \
    && Rscript -e 'devtools::install_github("JiekaiLab/dior")' \
               -e 'devtools::install_github("mojaveazure/seurat-disk")' \
               -e 'IRkernel::installspec()' \
    && rm -rf /opt/conda/pkgs/*
