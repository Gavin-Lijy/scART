\name{scART}
\alias{scART}
\title{scART for scATAC-seq}
\usage{
CreatescART(data,metadata = annotation)
RunImputation(art,k=1)
SparseFilter(art, ncell=2, ncell2=0.8, ncell3=2, nbin=10)
RunSim(art)
DimReduce(art)
set.seed(10)
RunTSNE(art, nSV=20, ndims=2, perplexity=30)
RunUMAP(art)
RunCluster(art,delta_cutoff = 4,rho_cutoff = 8)
RunTrajectory(art,anno='cell_type',nSV <- 20,ndim <- 10,gamma <- 10)
RunChromVAR(art,Org=c('hg19'),species = c("Homo sapiens")
MapBin2Gene(art,Org = 'hg19')
}
\description{
scART is a pipeline for single cell ATAC-seq (scATAC-seq).
}
\examples{
CreatescART(data,metadata = annotation)
}
