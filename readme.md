# rexseek

## Overview

R package for matrix processing in the exSeek project

## Install

```r
if (!requireNamespace("devtools")) install.packages('devtools');
devtools::install_github('dongzhuoer/rexseek');
```

## develop

Refer to this [post](https://dongzhuoer.github.io/_redirects/develop-upon-my-r-package.html)

## Non-Git

`data-raw/external/`

```r
'data-raw/external/all.txt' %>% readr::read_tsv() %>% 
	dplyr::select('transcript_id', gene_type = 'transcript_type') %>%
    readr::write_rds('data-raw/rna_type.rds', 'xz')
```

## Credits

`rna_type` comes from [Binbin Shi](https://github.com/ltbyshi)
