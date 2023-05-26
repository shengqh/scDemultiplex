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

# Speed (in seconds)

|                      | batch1_c1 | batch1_c2 | batch2_c1 | batch2_c2 | batch3_c1 | batch3_c2 |
| -------------------- | --------- | --------- | --------- | --------- | --------- | --------- |
| scDemultiplex_cutoff | 8.114     | 8.920     | 28.546    | 22.207    | 25.503    | 23.733    |
| scDemultiplex        | 227.507   | 159.049   | 249.854   | 248.873   | 394.811   | 682.445   |
| HTODemux             | 1.314     | 1.519     | 2.592     | 2.952     | 3.835     | 3.240     |
| MULTIseqDemux        | 0.660     | 0.589     | 1.104     | 1.203     | 1.475     | 1.458     |
| GMM_Demux            | 5.202     | 1.582     | 1.629     | 1.752     | 1.676     | 2.201     |
| BFF_raw              | 9.908     | 8.935     | 12.468    | 13.414    | 15.351    | 14.397    |
| BFF_cluster          | 17.532    | 18.368    | 33.113    | 35.795    | 46.347    | 47.997    |
| demuxmix             | 1.181     | 1.173     | 2.321     | 3.193     | 3.884     | 3.893     |
| hashedDrops          | 6.237     | 0.075     | 0.052     | 0.050     | 0.057     | 0.094     |

Both scDemultiplex_cutoff and scDemultiplex used mutli-thread. 
