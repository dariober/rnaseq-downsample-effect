<!-- vim-markdown-toc GFM -->

* [Set up](#set-up)
* [Usage](#usage)

<!-- vim-markdown-toc -->


Set up
======

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh

conda install mamba -n base -c conda-forge

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```


```
conda create --yes -n rnaseq-downsample-effect
conda activate rnaseq-downsample-effect
mamba install --yes --file requirements.txt -n rnaseq-downsample-effect
```


Usage
=====

```
snakemake -p -n -j 10 -C \
    sample_sheet=$PWD/sample_sheet.tsv \
    -d rnaseq-downsample-effect
```
