################################################
## pakeages ##
library(Matrix);library(proxy);library(gplots);library(Rtsne);
library(densityClust);library(scatterplot3d);library(Biobase);
library(RColorBrewer);library(ggplot2);
library(monocle);library(text2vec);
library(ChIPseeker)
library(GenomicFeatures)
library(data.table)
library(textTinyR)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg19)
library(R.matlab)
library(scales)

setwd("./Figures/scART_PACKAGES/0416_v2/")
source("./Figures/scART.R")
# source("/nethome/gjx/Project/scRNA-seq/scATAC-seq/GSE96769/HSC_scART_v2/Figures/function0416-gjx.R")
source("./Figures/function_0416_v2.R")

load("./ncounts.Rdata")
load("./annotation2_2347cells.Rdata")
dim(ncounts)
dim(annotation2)

annotation <- annotation2
art <- CreatescART(ncounts,metadata = annotation)  #metadata = annotation)
art@bmat$filter <- ncounts
head(art@metaData)
art <- RunSim(art)
art <- DimReduce(art, n=150, num=100, scale=F)

set.seed(10) 
art <- RunCluster(art, nSV=20, delta_cutoff = 11, rho_cutoff = 45)

set.seed(10) 
art <- RunTSNE(art, nSV=20, ndims=2, perplexity=30)
art <- RunUMAP(art, nSV=20)

tsnecols = c(brewer.pal(12, "Paired"), brewer.pal(7, "Accent"), brewer.pal(9, "Pastel1"))
tsnecols = unique(tsnecols)
pdf("Visualization_2D.pdf")
Visualization_2D(art,reductions = 'TSNE') 
Visualization_2D(art,reductions = 'UMAP')
Visualization_2D(art,reductions = 'TSNE',anno="type2",color=tsnecols,fileName="cell_type") 
Visualization_2D(art,reductions = 'UMAP',anno="type2",color=tsnecols) 
dev.off()

save(art, file="scART.Rdata")

load("scART.Rdata")


art <- RunTrajectory(art, nSV=20, anno="type2")

pdf("plotTrajectory.pdf")
plotTrajectory(art, anno="type2") 
dev.off()

art <- MapBin2Gene(art, ### the cell-by-bin matrix
                   binFormat = 'binary_matrix', ### the format of cell-by-bin matrix
                   bin_file = NULL,
                   Org = 'hg19', ### mm10
                   OrgDb = 'org.Hs.eg.db', ### 
                   TxDb = NULL, ### if Org = manual, you should input the TxDb defined by yourself 
                   convert_mat = TRUE, ### whether convert bin-by-cell matrix to cell-by
                   TSS_window = 5000 ### the window size around TSS to define the promoter 
                   )

p1 <- PlotSelectGenesATAC(art, gene2plot = c("GATA1"), reduction = 'TSNE', ncol = NULL)

pdf("PlotSelectGenesATAC.pdf")
p1
dev.off()



