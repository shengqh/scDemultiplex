scDemultiplex
==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Citation](#citation)
* [Usage](#example)
* [Speed test](#speed)
<a name="introduction"/>

# Introduction

Single-cell sequencing have been widely used to characterize cellular heterogeneity. Sample multiplexing where multiple samples are pooled together for single-cell experiments, attracts wide attention due to its benefits of increasing capacity, reducing costs, and minimizing batch effects. To analyze multiplexed data, the first crucial step is to demultiplex, the process of assigning cells to individual samples. Inaccurate demultiplexing will create false cell types and result in misleading characterization. We propose scDemultiplex, which models hashtag oligo (HTO) counts with beta-binomial distribution and uses an iterative strategy for further refinement. Compared with seven existing demultiplexing approaches, scDemultiplex achieves the highest accuracy in identifying singlets and removing negatives and multiplets.

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

Here is one [Example](http://htmlpreview.github.io/?https://github.com/shengqh/scDemultiplex/blob/main/vignettes/scDemultiplex.html).

For detailed analysis of the paper, please have a look at [scDemultiplex_analysis](https://github.com/shengqh/scDemultiplex_analysis).

# Speed test

| Dataset   | Computer              | CPU                                      | Memory | System   | R     | Cells | scDemultiplex cutoff | scDemultiplex refinement |
| --------- | --------------------- | ---------------------------------------- | ------ | -------- | ----- | ----- | -------------------- | -------------------- |
| batch1_c1 | ACCRE Cluster Gateway | Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz | 1000gb | CentOs 7 | 4.3.0 | 11900 | 8.1 sec              | 3.8 min              |
| batch1_c2 | ACCRE Cluster Gateway | Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz | 1000gb | CentOs 7 | 4.3.0 | 12923 | 8.9 sec              | 2.7 min              |
| batch2_c1 | ACCRE Cluster Gateway | Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz | 1000gb | CentOs 7 | 4.3.0 | 24905 | 28.5 sec             | 4.2 min              |
| batch2_c2 | ACCRE Cluster Gateway | Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz | 1000gb | CentOs 7 | 4.3.0 | 25763 | 22.2 sec             | 4.1 min              |
| batch3_c1 | ACCRE Cluster Gateway | Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz | 1000gb | CentOs 7 | 4.3.0 | 32886 | 25.5 sec             | 6.6 min              |
| batch3_c2 | ACCRE Cluster Gateway | Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz | 1000gb | CentOs 7 | 4.3.0 | 31956 | 23.7 sec             | 11.4 min             |

Both scDemultiplex cutoff mode and the following scDemultiplex refinement used mutli-thread. 
