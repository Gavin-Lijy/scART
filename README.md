# scART
![flow](image/flow.png)


# scART
**scART** is a R package for recognizing cell clusters and constructing trajectory from single-cell epigenomic data

## Introduction
scART is a Shen-lab in-house bioinformatics pipeline for single cell ATAC-seq (scATAC-seq).
scART is a statistical and powerfule toolkit to liminate noise and select significant features for cell state 
recognition and lineage construction from sparse single-cell epigenetic data. Using a compendium of single cell 
ATAC-seq datasets, scART exploited for robust recognizion of cell types and relevant regulons, as well as 
achieved learning developing trajectory from single-cell epigenetic data. scART decodes the regulatory heterogeneity 
in cell populations and reconstructed cell fate transition.

## Dependencies (for R >= 3.4.4)
The following packages have to be installed manually before installing scART:

The following packages have to be installed manually before installing scART:

```{r}
if (!requireNamespace(c("chromVAR","GenomicFeatures","GenomicRanges","motifmatchr","JASPAR2018","textTinyR","Matrix","text2vec","irlba","Rtsne","densityClust","scales","ggplot2","data.table","ChIPseeker","uwot","ggpubr","cowplot","SummarizedExperiment","monocle","RColorBrewer","scatterplot3d")),quietly = TRUE)
install.packages(c("chromVAR","GenomicFeatures","GenomicRanges","motifmatchr","JASPAR2018","textTinyR","Matrix","text2vec","irlba","Rtsne","densityClust","scales","ggplot2","data.table","ChIPseeker","uwot","ggpubr","cowplot","SummarizedExperiment","monocle","RColorBrewer","scatterplot3d"))
```

Now, you are now ready to install scART:

```{r}
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("JingxinGuo/scART") 
library(scART)
```

# Load tutorial data and annotation  
```{r}
setwd('/home/lijingyu/art/')
load('tutorial_data.Rdata')
```

# source the function
```{r}
source('scART.R')
source('function.R')
```

# Creating the scART object, replacing missing value and filtering bins
```{r message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE}
art <- CreatescART(data,metadata = annotation)  
# You can also use snap/10X output directory by Read_snap and Read_10X(Read_10X_h5) respectively.
# Or you can put the barcodes.tsv, bins.bed, matrix.mtx in a fold and provide Read_counts() with the fold address.

art <- RunImputation(art,k=1)
art <- SparseFilter(art, ncell=2, ncell2=0.8, ncell3=2, nbin=10)
```

# ![statistics](image/statistics.png)Dimensionality reduction 

```{r include=FALSE}
art <- RunSim(art)
art <- DimReduce(art)
```

![dim](README.assets/dim.png)

# Group cells into clusters

You can take a good look at the output pdf to adjust 'rho_cutoff' and 'delta_cutoff'

```{r message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE}
set.seed(10) 
art <- RunCluster(art,delta_cutoff = 4,rho_cutoff = 8)
```

# Visualize an Embedding
```{r}
set.seed(10) 
art <- RunTSNE(art, nSV=10, ndims=2, perplexity=30)
art <- RunUMAP(art)
p1 <- Visualization_2D(art,reductions = 'UMAP') 
p2 <- Visualization_2D(art,reductions = 'TSNE')
library(patchwork)
p1|p2

```

![umap](README.assets/umap.png)

# Trajectory analysis and visualization

```{r}
art <- RunTrajectory(art, anno='cell_type', nSV = 10, ndim= 3, gamma = 10)
plotTrajectory(art)
```

![trajectory](README.assets/trajectory.png)

# Create cell-by-gene matrix and explore gene accessibility score

```{r}
art <- MapBin2Gene(art, Org = 'hg19')
PlotSelectGenesATAC(art, gene2plot = c('GATA1','EBF1'))  
```

![plotGene](image/plotGene.png)

# RunChromVAR 


```{r}
art <- RunChromVAR(art,Org=c('hg19'),species = c("Homo sapiens") ,min.count=10)
PlotSelectTF(art,TF2plot=c('GATA1'))
### min.count:the threshold of a peaks found at at least 10 cells
```

![motif](image/motif.png)

```R
save(art,'scart.Rdata')
```
