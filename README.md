
<!-- README.md is generated from README.Rmd. Please edit that file -->

# INDEED

## Overview

This package implements INDEED algorithm from Zuo et. al.’s Methods
paper: INDEED: Integrated differential expression and differential
network analysis of omic data for biomarker discovery (PMID: 27592383).

This R package will generate a csv file containing information such as
p-values, node degree and activity score for each biomolecule. A higher
activity score indicates that the corresponding biomolecule has more
neighbors connceted in the differential network and their p-values are
more statistically significant. It will also generate a csv file for the
differential network created by INDEED.

## Installation

You can install INDEED from github with:

``` r
# The development version from GitHub:
# install.packages("devtools")
devtools::install_github("cx30/INDEED")
```

## Usage

``` r
library(INDEED)

# Example 1:
# Using partial correlation to obtain sparse differential network
select_sig(Met_GU, Met_Group_GU, Met_name_GU, partial = TRUE)

# Example 2:
# Using Spearman correlation to obtain differential network
select_sig(Met_GU, Met_Group_GU, Met_name_GU, partial = FALSE, method = "spearman")
```
