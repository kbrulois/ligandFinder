## Overview

ligandFinder is a tool to run large-scale screening to identify GPCR ligands

## Installation

From within R:

```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")

if (!requireNamespace("arrow", quietly = TRUE)) install.packages("arrow")
install.packages("bio3d", dependencies=TRUE)
install.packages("httr", dependencies=TRUE)

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("kbrulois/ligandFinder")
```

## Tutorial

Please see package website for full documentation:

https://kbrulois.github.io/ligandFinder....

