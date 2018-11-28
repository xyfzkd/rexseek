# rexseek

## Overview

R package for matrix processing in the exSeek project

## Install

Some functions only works in Bioconductor 3.7

> For Students in Lu Lab, this step is not needed (and even causes error) if you run R on cnode.

```r
if (!requireNamespace("BiocManager")) install.packages('BiocManager');
BiocManager::install(version = '3.7', ask = F, update = F)
```

Then you can install the package

```r
if (!requireNamespace("remotes")) install.packages('remotes');
remotes::install_github('dongzhuoer/rexseek', upgrade = F);
```

## Usage

For students in Lu Lab, please refer to [this guide](https://github.com/dongzhuoer/lulab-rotation-summary/blob/master/exseek.md)

## develop

Refer to this [post](https://dongzhuoer.github.io/_redirects/develop-upon-my-r-package.html)

## Non-Git


### files neccessary for developing the packages

```r
'data-raw/external/all.txt' %>% readr::read_tsv() %>% 
	dplyr::select('transcript_id', gene_type = 'transcript_type') %>%
    readr::write_rds('data-raw/rna_type.rds', 'xz')
```

### files used in real-case testing

`data-raw/external/` can be found in `dongzhuoer/lulab-rotation-summary` repo, `exseek/` folder

## Credits

`rna_type` comes from [Binbin Shi](https://github.com/ltbyshi)
