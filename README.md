# 2022 DIGIT-BIO workshop: "Do my clusters make sense? Some statistical and biological clues to help..."

<img src="digitbio_logo.png" align="center" />

Unsupervised classification (clustering) is used in many [fields](https://digitbio-ia.github.io/sequences/concepts/s2_clustering). But once a clustering result has been obtained, how can we evaluate its pertinence? In this work, we will work with several statistical and biological indices to reinforce our confidence (or not...) in a clustering result. We will use as an example the identification of co-expressed genes from transcriptomic data with finite mixture models. 

This repo contains all materials for the tutorial. It is organized as follows:

- `biblio/`: folder containing a few helpful and relevant references
- `code/`: folder containing R code contained in the Rmarkdown tutorial
- `data/`: folder containing three `.rds` files that are to be downloaded and loaded in the tutorial R code. These files were generated with the code indicated in the Appendix of the Rmarkdown tutorial
- `doc/`: compiled Rmarkdown tutorial for rendering on Github pages 
- `slides/`: folder containing slides presented during the tutorial.

The rendered Rmarkdown HTML document can be viewed [here](https://www.andrea-rau.com/2022_DIGIT-BIO_workshop/).

To install the packages needed to run the tutorial, execute the following R code (~15 min of installation time):

```
install.packages("BiocManager")
BiocManager::install(c("mclust", "clValid", "fpc", "ggplot2",
                       "coseq", "clusterProfiler", "org.Mm.eg.db"))
```

In some cases, installation of the `clusterProfiler` package may cause
errors due to a dependence on the `r Biocpkg("ggtree")` package (see [here](https://github.com/YuLab-SMU/ggtree/issues/544). 
In this case, you can try
using the following code to install the development version of `ggtree`
before installing the remaining packages as described above.

```
install.packages("remotes")
remotes::install_github('YuLab-SMU/ggtree')
```