# rexseek

## Overview

R package for matrix processing in the exSeek project

## Install

```r
if (!requireNamespace("remotes")) install.packages('remotes');
remotes::install_github('dongzhuoer/rexseek');
```

## Usage

For Students in Lu Lab, please refer to [this guide](https://github.com/dongzhuoer/lulab-rotation-summary/blob/master/exseek.md)

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
