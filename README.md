scDemultiplex
==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Citation](#citation)
* [Usage](#example)
* [Speed test](#speed)
<a name="introduction"/>

# Introduction

scDemultiplex is an R package for ...

<a name="installation"/>

# Installation


```R
library(devtools)
install_github("shengqh/deMultiplex")
```

or

```R
BiocManager::install("shengqh/deMultiplex")
```
  
## Trouble shooting


<a name="citation"/>

# Citation

scDemultiplex: 

<a name="example"/>

# Usage

After installing scDemultiplex, use following codes to run examples:

## Example 1:


<a name="speed"/>

# Speed test

|Computer|CPU|Memory|System|R version|Example1|Example2|Example3|
|---|---|---|---|---|---|---|---|
|MacBook Pro Laptop|1.7 GHz Intel Core i5|4 Gb|MacOS|3.4.3|39 min |13 min|3 min|
|MacBook Pro Laptop|2.7 GHz Intel Core i5|8 Gb|MacOS|3.4.1|6.5 min|12.5 min|1.1 min|
|MacBook Pro Laptop|3.1 GHz Intel Core i5|16 Gb|MacOS|3.6.0|2.5 min|3.6 min|41 sec|
|Lenovo Laptop|1.8 GHz Intel(R) i7 8565U|16 Gb|Windows 10|3.6.0|3.4 min|5.1 min| 1.1 min|
|Windows Desktop|2 GHz Intel(R) Xeon(R) E5-2620|32 Gb|Windows 7|3.4.3|3 min|5 min| 1 min|
|Windows Desktop|2.6 GHz Intel(R) Xeon(R) E5-2640|64 Gb|Windows 10|3.6.0|4.1 min|6 min| 1.3 min|
|ACCRE Cluster Gateway|Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz|1000 Gb|CentOs 7|3.5.1|3.1 min|5.2 min| 28.6 sec|

