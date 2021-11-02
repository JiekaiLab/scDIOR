FROM hjfeng/base-jupyter-notebook:1.1
ARG DEBIAN_FRONTEND=noninteractive # 用于把参数设置默认为True，不用交互提供指令
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
        r-leidenbase \
        r-monocle3 \
    && pip install \
        scanpy==1.7.2 \
        scvelo==0.2.3 \
        diopy \
        loompy \
    && Rscript -e 'devtools::install_github("JiekaiLab/dior")' \
    && Rscript -e 'IRkernel::installspec()' \
    && rm -rf /opt/conda/pkgs/*
