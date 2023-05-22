scDemultiplex
==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Citation](#citation)
* [Usage](#example)
* [Speed test](#speed)
<a name="introduction"/>

# Introduction

Single-cell sequencing have been widely used to characterize cellular heterogeneity. Sample multiplexing where multiple samples are pooled together for single-cell experiments, attracts wide attention due to its benefits of increasing capacity, reducing costs, and minimizing batch effects. To analyze multiplexed data, the first crucial step is to demultiplex, the process of assigning cells to individual samples. Inaccurate demultiplexing will create false cell types and result in misleading characterization. We propose scDemultiplex, which models hashtag oligo (HTO) counts with beta-binomial distribution and uses an iterative strategy for further refinement. Compared with five existing demultiplexing approaches, scDemultiplex achieves the highest accuracy in identifying singlets and removing negatives and multiplets.

<a name="installation"/>

# Installation

```R
library(devtools)
install_github("shengqh/cutoff")
install_github("shengqh/scDemultiplex")
```

or

```R
BiocManager::install("shengqh/cutoff")
BiocManager::install("shengqh/scDemultiplex")
```

# Usage

Here is one [Example](vignettes/scDemultiplex.html).

For detailed analysis of the paper, please have a look at [scDemultiplex_analysis](https://github.com/shengqh/scDemultiplex_analysis).
