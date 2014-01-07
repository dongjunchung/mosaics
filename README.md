mosaics
=======

MOSAiCS (MOdel-based one and two Sample Analysis and Inference for ChIP-Seq) is a flexible statistical framework for  the one- (ChIP sample) and two-sample (ChIP sample and matched control sample) analysis of ChIP-seq data. In addition, MOSAiCS-HMM extends MOSAiCS with a Hidden Markov Model (HMM) structure to account for spatial dependence and call broad peaks as in the case of histone modifications. mosaics package provides computationally efficient and user friendly interface to process ChIP-seq data, implement exploratory analysis, fit MOSAiCS and MOSAiCS-HMM models, call peaks, and export peak lists for downstream analysis.

Stable versions of mosaics package is maintained through Bioconductor. To install or update the stable version of mosaics package, please run:

```
source("http://bioconductor.org/biocLite.R")
biocLite("mosaics")
```

MOSAiCS vignette provides a good start point for the ChIP-seq data analysis using mosaics package and it can be found at http://www.bioconductor.org/packages/2.13/bioc/vignettes/mosaics/inst/doc/mosaics-example.pdf. Please check http://groups.google.com/group/mosaics_user_group for discussions and questions regarding ChIP-seq data analysis using mosaics package. You can track development of mosaics package at http://github.com/dongjunchung/mosaics.

Development
===========

To install the development version of mosaics, it's easiest to use the devtools package:

```
#install.packages("devtools")
library(devtools)
install_github("dongjunchung/mosaics")
```

References
==========

Kuan PF, Chung D, Pan G, Thomson JA, Stewart R, and Keles S (2011), "A statistical framework for the analysis of ChIP-Seq data," _Journal of the American Statistical Association_, 106: 891-903.

Chung D, Zhang Q, Keles S (2013), "MOSAiCS-HMM: A model-based approach for detecting regions of histone modifications from ChIP-seq data," Datta S and Nettleton D (eds.), _Statistical Analysis of Next Generation Sequence Data_, Springer.
