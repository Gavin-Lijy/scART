   # library(methods)
  # methods::setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
  # setClass('scART',slots=list(barcode="character",feature='GRanges',metaData="data.frame",bmat = "list",smat='Matrix',gmat = "Matrix",
  #                             mmat = "Matrix",reductions = "list",trajectory='MatrixOrmatrix' ))
  # 
  # .valid.scART.barcode <- function(object)
  # {
  #   if(length(object@barcode) != nrow(object@metaData)){
  #     return("slot 'barcode' have different length from 'metaData'")		
  #   }
#   NULL
# }
# 
# .valid.scART <- function(object)
# {
#   
#   c(.valid.scART.barcode(object))
# }
# # methods::setValidity("scART", .valid.scART)

setClass( 'SVD',slots = c(x='matrix',sdev='numeric'))
setClass('TSNE',slots = c(matrix='matrix',nSV='numeric'))
setClass('TSNE_3D',slots = c(matrix='matrix',nSV='numeric'))


CreatescART=function(data,barcode,bins,metadata){
  if(missing(data)){
    stop("data is missing");
  }
  if(missing(barcode)||missing(bins)){
    barcode <- colnames(data)
    bins <- rownames(data)
  }
  library(textTinyR)
  barcode<-as.character(barcode)
  obj<-new('scART')
  if(missing(metadata)){
    depth <- sparse_Sums(data, rowSums = F)
    obj@metaData<-data.frame(depth)
    rownames(obj@metaData)<-barcode
  }else{obj@metaData<-metadata}
  
  obj@bmat=list(NULL,NULL,NULL,NULL,NULL)
  names(obj@bmat)=c('raw','binary','imputation','filter','TF_IDF')
  obj@bmat$raw=data
  obj@barcode<-barcode
  # bins<-lapply(bins, function(x) {strsplit(x,split = '-')})
  bins<-rownames(data)
  bins_bed = do.call(rbind, strsplit(x = bins, split = '[:-]'))
  df<- data.frame(bins_bed,stringsAsFactors=FALSE)
  colnames(df)<-c('chr','start','end') 
  library(GenomicRanges)
  features<-makeGRangesFromDataFrame(df)
  obj@feature<-features
  binarize= function(x,threshold=NA) {
    
    x@x[x@x <= threshold] = 0
    x@x[x@x > threshold] = 1
    return(x)
  }
  obj@bmat$binary<-binarize(data,threshold = 0)
  nCounts<-sparse_Sums(data, rowSums = F)
  obj@metaData$nCounts<-nCounts
  return(obj)
}

RunImputation=function(obj,k=1,ratio=1){
  library(Matrix)
  m=obj@bmat$raw
  knn=function(m,k=1,ratio=1){
    binarize= function(x,threshold=NA) {
      x@x[x@x < threshold] = 0
      x@x[x@x >=threshold] = 1
      return(x)
    }
    rows<-nrow(m)
    if(!k%in%1:5){stop('k should be int in 1:5')}
    new<-m[1:(rows-2*k),]+m[(2*k+1):rows,]
    rownames(new)<-rownames(m)[(k+1):(rows-k)]
    if(k>1){ 
      for (i in (k-1):1 ){ 
        new<-new+m[(1+i):(rows-2*k+i),]+m[(2*k+1-i):(rows-i),]
      }}
    
    bi<-binarize(new,threshold=2*k*ratio)
    m_sum<-bi+m[(k+1):(rows-k),]
    matrix<-binarize(m_sum,threshold=0.5)
    final<-rbind(m[1:k,],matrix,m[(rows-k+1):rows,])   
    final@Dimnames[[1]]=m@Dimnames[[1]]
    final<-as(as.matrix(final),"Matrix")
    return(final)
    
  }
  final<-knn(m,k=k,ratio = ratio)
  
  obj@bmat$imputation<-final
  return(obj)
}


SparseFilter <- function(obj, ncell=NULL, ncell2=NULL, ncell3=NULL,nbin=NULL, genome='hg19') {
  raw <- obj@bmat$imputation
  if(genome=='hg19'){
    system('ls')
    system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/seq.cov01.ONHG19.bed.gz");
    library(GenomicRanges);
    black_list = read.table("seq.cov01.ONHG19.bed.gz");}
  if(genome=='mm10'){system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz");
    library(GenomicRanges);
    black_list = read.table("mm10.blacklist.bed.gz");}
    black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
  );
  idy = queryHits(findOverlaps(obj@feature, black_list.gr));
  if(length(idy) > 0){raw <- raw[-idy,]}
  
  
  if (is.null(ncell)) {
    ncell <- dim(raw)[2]*0.05
  } else {
    ncell <- ncell
  }
  if (is.null(ncell2)) {
    ncell2 <- dim(raw)[2]*0.75
  } else {
    ncell2 <- dim(raw)[2]*ncell2
  }
  if (is.null(ncell3)) {
    ncell3 <- 2
  } else {
    ncell3 <- ncell3
  }
  
  npeak=nbin
  if (is.null(npeak)) {
    npeak <- 50
  } else {
    npeak <- npeak
  }
  ncounts <- raw
  library(textTinyR)
  new_peaks = sparse_Sums(raw, rowSums = T)
  new_counts = sparse_Sums(raw, rowSums = F)
  ncounts <- raw[new_peaks > ncell & new_peaks <= ncell2,] 
  new_peaks2 = sparse_Sums(ncounts, rowSums = T)
  new_peaks2 <- scale(log10(new_peaks2))
  pdf("Cell_Site_Distribution.pdf", width=10)    
  par(mfrow=c(1,2))
  options(repr.plot.width=4, repr.plot.height=4)
  hist(log10(new_peaks+1),main="No. of Cells Each Site is Observed In",breaks=50,xlab="log10, nCells")
  # abline(v=log10(ncell),lwd=2,col="indianred")
  # abline(v=log10(ncell2),lwd=2,col="indianred")
  hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50,xlab="log10, nBins")
  # abline(v=log10(npeak),lwd=2,col="indianred")
  # hist(scale(log10(new_peaks+1)),main="No. of Cells Each Site is Observed In",breaks=50,xlab="raw, scale(log10())")
  # abline(v=ncell3,lwd=2,col="indianred")
  # hist(scale(log10(new_peaks[new_peaks2 < ncell3]+1)),main="No. of Cells Each Site is Observed In",xlab="filter1,scale(log10())",breaks=50)
  # abline(v=ncell3,lwd=2,col="indianred")
  dev.off()
  
  ncounts2 <- ncounts[which(new_peaks2 < ncell3),]
  cell.use <- colnames(ncounts2)
  ncounts2 <- ncounts2[,new_counts >= npeak]
  ncounts <- ncounts2
  new_peaks = sparse_Sums(ncounts, rowSums = T)
  # pdf("Cell_Site_Distribution_after_filter.pdf", width=10)    
  # par(mfrow=c(1,2))
  # options(repr.plot.width=4, repr.plot.height=4)
  # hist(log10(new_peaks+1),main="No. of Cells Each Site is Observed In",breaks=50,xlab="log10()")
  # hist(scale(log10(new_peaks+1)),main="No. of Cells Each Site is Observed In",breaks=50,xlab="scale(log10())")
  # dev.off()
  obj@bmat$filter<-ncounts
  
  colname<-colnames(obj@metaData)
  obj@metaData<-as.data.frame(obj@metaData[colnames(ncounts),])
  rownames(obj@metaData)<-colnames(ncounts)
  colnames(obj@metaData)<-colname
  obj@barcode<-colnames(ncounts)
  bins <- rownames(obj@bmat$filter)
  bins_bed = do.call(rbind, strsplit(x = bins, split = '[:-]'))
  bins_bed = as.data.frame(bins_bed)
  write.table(bins_bed, "./bins_bed.bed", quote=F, sep="\t",row.names=F, col.names=F)
  peakfile <- "./bins_bed.bed"
  bin.use <- chromVAR::getPeaks(peakfile, sort_peaks = TRUE)
  file.remove("./bins_bed.bed")
  idy = queryHits(findOverlaps(art@feature, bin.use));
  features <- art@feature
  if(length(idy) > 0){features <- features[-idy,]}  
  obj@feature<-features
  return(obj)
}

RunSim <- function(obj) {
  library(textTinyR)
  
  ncounts <- obj@bmat$filter
  col <- sparse_Sums(ncounts, rowSums = F)
  summary(col)
  row <- sparse_Sums(ncounts, rowSums = T)
  summary(row)
  a <- t(as.matrix(ncounts) )
  nfreqs = t(a / col)
  idf = as(log(1 + ncol(ncounts) / row), "sparseVector")
  a <- as(Diagonal(x=as.vector(idf)), "sparseMatrix")
  tf_idf_counts = a %*% nfreqs
  obj@bmat$TF_IDF<-tf_idf_counts
  transformed <- tf_idf_counts
  library(text2vec)
  data <- t(transformed)
  cosine_sim <- sim2(data ,data , method = "cosine",norm="l2")
  obj@smat<-cosine_sim
  return(obj)
}

DimReduce <- function(obj, n=NULL,span=NULL,num=NULL,scale=NULL) {
  library(irlba)
  data<-obj@smat
  if (is.null(data)) {
    data<-cellSim
  } else {
    data<-data
  }
  
  if (is.null(n)) {
    n<-50
  } else {
    n<-n
  }
  if (is.null(span)) {
    span<-0.5
  } else {
    span<-span
  }
  if (is.null(num)) {
    num<- 50
  } else {
    num<-num
  }
  if (is.null(scale)) {
    scale <- T
  } else {
    scale <- scale
  }
  
  n_svd=n
  set.seed(0)
  svd = prcomp_irlba(data, n = n_svd, scale= scale, center=T)
  decomp_data  <- svd$x 
  svd.var <- svd$sdev
  
  mydf <- data.frame(x=1:length(svd$sdev[c(1:num)]), y=svd$sdev[c(1:num)])
  ret <- loess(y~x, data=mydf, span=span)
  newX=seq(1,length(svd$sdev[c(1:num)]),1)
  ret_p <- data.frame(x=newX, y=predict(ret, newdata=data.frame(x=newX)))
  pdf("Standard Deviation of SVs.pdf")
  x =c(1:num)
  plot(x,svd$sdev[1:num], pch=16,xlab="SVs",ylab="Standard Deviation of SVs")
  # lines(newX[1:num],ret_p$y[1:num], pch=16,col="red")
  dev.off()
  
  svd_obj<-new('SVD')
  svd_obj@x<-as.matrix(svd$x )
  svd_obj@sdev<-svd$sdev
  obj@reductions=list(NULL,NULL,NULL,NULL)
  names(obj@reductions)=c('SVD','TSNE','TSNE_3D','UMAP')
  obj@reductions$SVD<-svd_obj
  
  return(obj)
}


RunTSNE <- function(obj, nSV=NULL, ndims=NULL, perplexity =NULL,seed.use=NULL) {
  
  library(Rtsne)
  if (is.null(ndims)) {
    ndims<-2
  } else {
    ndims<-ndims
  }
  if (is.null(nSV)) {
    nSV<-15
  } else {
    nSV <- nSV
  }
  if (is.null(perplexity)) {
    perplexity<-30
  } else {
    perplexity <- perplexity
  }
  if (is.null(seed.use)) {
    seed.use <- 10
  } else {
    seed.use <- seed.use
  }
  
  decomp_data  <- obj@reductions$SVD@x 
  pca.var <- obj@reductions$SVD@sdev
  set.seed(seed.use)
  
  cut <- pca.var[1]/pca.var[2]
  if (cut > 2){
    svd_tsne = decomp_data[,2:nSV]
  }else{
    svd_tsne = decomp_data[,1:nSV]
  }
  print(dim(svd_tsne)[2])
  
  set.seed(10)
  data <- svd_tsne
  tsne_out = Rtsne(data, pca=F, dims = ndims, perplexity = perplexity,
                   theta = 0.5, check_duplicates = F, max_iter = 500)
  tsne_out <- tsne_out$Y
  if (ndims==2) {
    colnames(tsne_out) <- c("tSNE_1","tSNE_2")
  } else {
    colnames(tsne_out) <- c("tSNE_1","tSNE_2","tSNE_3")
  } 
  
  rownames(tsne_out) <- obj@barcode
  if(ndims==2){
    TSNE<-new('TSNE')
    TSNE@matrix<-tsne_out
    TSNE@nSV<-nSV
    obj@reductions$TSNE<-TSNE
  }else{TSNE_3D<-new('TSNE_3D')
  TSNE_3D@matrix<-tsne_out
  TSNE_3D@nSV<-nSV
  obj@reductions$TSNE_3D<-TSNE_3D}
  return(obj) 
  
}

RunCluster <- function(obj,rho_cutoff,delta_cutoff,tsne_3D,nSV) {
  
  if(missing(rho_cutoff)){
    rho_cutoff=2
  }
  if(missing(delta_cutoff)){
    delta_cutoff=10
  }
  library(Rtsne)
  if(missing(nSV)){
    nSV<-15}
  
  if (missing(tsne_3D)) {
    
    decomp_data  <- obj@reductions$SVD@x 
    pca.var <- obj@reductions$SVD@sdev
    
    cut <- pca.var[1]/pca.var[2]
    if (cut > 2){
      svd_tsne = decomp_data[,2:nSV]
    }else{
      svd_tsne = decomp_data[,1:nSV]
    }
    print(dim(svd_tsne)[2])
    
    set.seed(10)
    
    data <- svd_tsne
    tsne_out = Rtsne(data, pca=F, dims = 3, perplexity = 30,
                     theta = 0.5, check_duplicates = TRUE, max_iter = 500)
    tsne_3D <- tsne_out$Y} else {
      tsne_3D <- tsne_3D
    }
  library(densityClust)
  
  set.seed(10)
  points <- tsne_3D
  dis <- dist(points)
  dclust <- densityClust(dis, gaussian=TRUE)
  dclust <- findClusters(dclust, 
                         rho=rho_cutoff, 
                         delta= delta_cutoff)
  
  pdf("Dclust_rho_delt.pdf")
  options(repr.plot.width=6, repr.plot.height=6)
  plot(dclust$rho,dclust$delta,pch=20,cex=0.6)
  points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20,cex=0.8)
  text(dclust$rho[dclust$peaks]+0.5,dclust$delta[dclust$peaks]+0.5,labels=dclust$clusters[dclust$peaks])
  abline(v=rho_cutoff)
  abline(h=delta_cutoff)
  dev.off()
  # save(dclust, file=paste(out_dir, "dclust.Rdata",sep=""))
  print(paste("number of clusters=",length(table(dclust$clusters))))
  scART_cluster <- as.factor(dclust$clusters)
  names(scART_cluster) <- obj@barcode
  
  obj@metaData$cluster<-scART_cluster
  return(obj)
}


Visualization_2D <- function(obj,reductions='TSNE' ,anno=NULL,fileName=NULL,color=NULL,size=NULL){
  library( scales)
  if (is.null(anno)) {
    anno <- obj@metaData$cluster
  } else {
    anno <- eval(parse(text =paste0('obj@metaData$',anno) ))
  }
  if (is.null(fileName)) {
    fileName <- "Cluster"
  } else {
    fileName <- fileName
  }
  if (is.null(color)) {
    tsnecols = hue_pal()(length(table(anno)))
    tsnecols = unique(tsnecols)
  } else {
    tsnecols = color
  }
  if (is.null(size)) {
    size <- 1
  } else {
    size<- size
  }
  if(reductions=='UMAP'){data<-as(obj@reductions$UMAP,'matrix');x = "UMAP_1"; y = "UMAP_2"}else{
    data<-obj@reductions$TSNE@matrix;x = "tSNE 1";y = "tSNE 2"}
  clusters <- anno
  points <- as.data.frame(data)
  # pdf(paste(fileName,"2D_Visualization.pdf",sep="_"),width=6,height=6)    
  data_df <- as.data.frame(points)
  sample_state <- clusters
  library(ggplot2)
  p <- ggplot(data = data_df, aes(x = data_df[,1], y =  data_df[,2]))
  p +  geom_point(aes(color = as.factor(sample_state)), size = 0.5, position=position_jitter(width=0.5,height=.25)) + 
    scale_color_manual(values=tsnecols) +
    labs(x = x, y = y, title="2D visualization",
         colour = fileName) +
    theme(legend.position = "right", legend.key.height = grid::unit(0.35, "in")) + 
    theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white", colour = "black"))
  
}


# Visualization_3D <- function(data, anno,fileName=NULL){
#   if (is.null(anno)) {
#     anno <- obj@metaData$cluster
#   } else {
#     anno <- anno
#   }
#   if (is.null(fileName)) {
#     fileName <- "Cluster"
#   } else {
#     fileName <- fileName
#   }
#   
#   clusters <- anno
#   points <- data
#   pdf(paste(fileName,"3D_Visualization.pdf",sep="_"),width=6,height=6)    
#   mydata<-cbind(points[,1],points[,2],points[,3])
#   colnames(mydata)<-c("1m","2m","3m")
#   mydata<-as.data.frame(mydata)
#   tsnecols = c(brewer.pal(12, "Paired"), brewer.pal(7, "Accent"), brewer.pal(9, "Pastel1"))
#   tsnecols = unique(tsnecols)
#   library(scatterplot3d)
#   with(mydata, {
#     scatterplot3d(points[,1],   # x axis
#                   points[,2],     # y axis
#                   points[,3], 
#                   color=tsnecols[as.factor(clusters)],
#                   main="3D visualization",
#                   pch=19,cex.symbols=0.4,xlab="tSNE 1",ylab="tSNE 2",zlab="tSNE 3",cex.axis=1)
#   })
#   legend("topright",pch=c(rep(20,length(names(table(clusters))))),col=tsnecols[1:length(table(as.factor(clusters)))],
#          legend=names(table(clusters)))
#   dev.off()
# }
# 
# 


# ### ******* 1. a function to map ATAC bins/peaks to genome annotation (gene levels) --------
# MapBin2Gene = function(Bmat = NULL, ### the cell-by-bin matrix
#                        binFormat = 'binary_matrix', ### the format of cell-by-bin matrix
#                        bin_file = NULL,
#                        Org = 'mm10', ### mm10
#                        OrgDb = 'org.Mm.eg.db', ### 
#                        TxDb = NULL, ### if Org = manual, you should input the TxDb defined by yourself 
#                        convert_mat = TRUE, ### whether convert bin-by-cell matrix to cell-by
#                        TSS_window = 3000 ### the window size around TSS to define the promoter 
# ){Bmat<-obj@bmat

# options(stringsAsFactors = F)

# cat(">> checking depedent packages ... \t\t\t", format(Sys.time(), 
#                                                        "%Y-%m-%d %X"), "\n")

# if (!(requireNamespace("BiocManager", quietly = TRUE))) {
#   cat("Please install package 'BiocManager'");
#   install.packages('BiocManager')
# }
# library(data.table)
# library(GenomicFeatures)
# library(ChIPseeker)

# if (!(requireNamespace("data.table", quietly = TRUE))) {
#   cat("Please install package 'data.table'");
#   install.packages('data.table')
# }

# if (!(requireNamespace("GenomicFeatures", quietly = TRUE))) {
#   cat("Please install package 'GenomicFeatures'");
#   BiocManager::install('GenomicFeatures')
# }

# if (!(requireNamespace("ChIPseeker", quietly = TRUE))) {
#   cat("Please install package 'ChIPseeker'");
#   BiocManager::install('ChIPseeker')
# }



# ### load the ATAC IDs of peaks or bins -------  
# cat(">> loading ATAC features ...\t\t\t", format(Sys.time(), 
#                                                  "%Y-%m-%d %X"), "\n")

# if(sum(!is.matrix(Bmat))>0){
#   Bmat = Matrix::as.matrix(Bmat)
#   ### matrix is faster than data.frame
# }

# if(binFormat == 'binary_matrix'){
#   bins = as.character(rownames(Bmat))
#   bins_bed = do.call(rbind, strsplit(x = bins, split = '[:-]'))
#   bins_bed = as.data.frame(bins_bed);colnames(bins_bed)[1:3] = c('chr','start','end')
#   bins_bed[,2] = as.numeric(bins_bed[,2])
#   bins_bed[,3] = as.numeric(bins_bed[,3])
#   bins.gr = GenomicRanges::GRanges(bins_bed[,1], IRanges(start = bins_bed[,2], end = bins_bed[,3], names = bins))
# }
# if(binFormat == 'bin_file'){
#   bins_bed = read.table(bin_file, sep='\t', fill =T, stringsAsFactors = F, header = F)
#   colnames(bins_bed)[1:3] = c('chr','start','end')
#   bins_bed = transform(bins_bed, name = paste0(chr,':',start,'-',end))
#   bins.gr = GenomicRanges::GRanges(bins_bed$chr, IRanges(start = bins_bed$start, end = bins_bed$end, names = bins))
# }


# ### load the genome reference  -------- 
# cat(">> adding gene annotation...\t\t\t", format(Sys.time(), 
#                                                  "%Y-%m-%d %X"), "\n")

# if(Org == 'mm9'){
#   if (!(requireNamespace("TxDb.Mmusculus.UCSC.mm9.knownGene", quietly = TRUE))) {
#     message("Please install package 'TxDb.Mmusculus.UCSC.mm9.knownGene'");
#     BiocManager::install('TxDb.Mmusculus.UCSC.mm9.knownGene')
#   }
#   txdb = TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
#   OrgDb = 'org.Mm.eg.db'
# }
# if(Org == 'mm10'){
#   if (!(requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE))) {
#     message("Please install package 'TxDb.Mmusculus.UCSC.mm10.knownGene'");
#     BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene')
#   }
#   txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
#   OrgDb = 'org.Mm.eg.db'
# }
# if(Org == 'hg19'){
#   if (!(requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE))) {
#     message("Please install package 'TxDb.Hsapiens.UCSC.hg19.knownGene'");
#     BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')
#   }
#   +{}
#   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#   OrgDb = 'org.Hs.eg.db'
# }
# if(Org == 'hg38'){
#   if (!(requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))) {
#     message("Please install package 'TxDb.Hsapiens.UCSC.hg38.knownGene'");
#     BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')
#   }
#   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#   OrgDb = 'org.Hs.eg.db'
# }
# if(Org == 'manual'){
#   txdb = TxDb
# }

# ### map to genome reference -----
# cat(">> mapping bins to genome...\t\t\t", format(Sys.time(), 
#                                                  "%Y-%m-%d %X"), "\n")
# TSS_window <- TSS_window
# binAnno <- ChIPseeker::annotatePeak(bins.gr, tssRegion=c(-TSS_window, TSS_window),
#                                     TxDb=txdb, annoDb=OrgDb, level = 'gene',
#                                     genomicAnnotationPriority = c("Promoter","Exon", "Intron",
#                                                                   "Downstream", "Intergenic")) 

# bin2gene = as.data.frame(binAnno@anno)
# setDT(bin2gene,keep.rownames = 'ID')
# head(bin2gene)

# bin2gene[,Location:=annotation]
# bin2gene[grep('Exon',annotation)]$Location = 'Exon'
# bin2gene[grep('Intron',annotation)]$Location = 'Intron'
# bin2gene[grep('Promoter',annotation)]$Location = 'Promoter'
# bin2gene[grep('Down',annotation)]$Location = 'Downstream'

# bin2use = bin2gene[grep("Promoter",Location)]
# gene2use = unique(bin2use$SYMBOL)

# ### select the ATAC bins in gene's promoter
# if(convert_mat==TRUE){
#   ### convert cell-by-bin matrix to cell-by-gene matrix 
#   mat2plot = do.call(rbind, lapply(gene2use, function(gene){
#     bins2use = bin2use[SYMBOL==gene]$ID
#     if(length(bins2use)==1){
#       sum_bins = Bmat[bins2use,]
#     }else{
#       sum_bins = colSums(Bmat[bins2use,])
#     }
#     sum_bins
#   }))
#   rownames(mat2plot) = gene2use
#   cat(">>", nrow(mat2plot), "genes with promoter ATAC bins retained after filtering ...\t\t\t", format(Sys.time(), 
#                                                                                                        "%Y-%m-%d %X"), "\n") 
#   gmat <- as(mat2plot,'Matrix')
#   # gmat <- t(t(gmat)/art@metaData$nCounts)
#   obj@gmat<-gmat
#   return(obj)
# }else{
#   cat(">>", length(gene2use), "genes with promoter ATAC bins retained after filtering ...\t\t\t\n", 'finished in ',format(Sys.time(), 
#                                                                                                                           "%Y-%m-%d %X"), "\n")  
#   return(bin2gene)
# }
# }

#### ***** Code examples:
### using the downloaded genome TxDb:
### gene_matrix = MapBin2Gene(Bmat = binary_matrix,Org = 'mm10', convert_mat = TRUE, TSS_window = 3000)

### using manually constructed TxDb:
### TxDb = makeTxDbFromGFF(gtf_file);
### gene_matrix = MapBin2Gene(Bmat = binary_matrix,Org = 'manual', OrgDb = 'org.Mm.eg.db', TxDb = TxDb, convert_mat = TRUE, TSS_window = 3000)
RunUMAP<-function(obj,dims=2,seed.use=10,nSV=20){
  library(uwot)
  decomp_data  <- obj@reductions$SVD@x 
  pca.var <- obj@reductions$SVD@sdev
  set.seed(seed.use)
  
  cut <- pca.var[1]/pca.var[2]
  if (cut > 2){
    svd_tsne = decomp_data[,2:nSV]
  }else{
    svd_tsne = decomp_data[,1:nSV]
  }
  print(dim(svd_tsne)[2])
  
  
  
  data <- svd_tsne
  umaps  <- uwot::umap(data,n_components = dims) 
  if(dims==2){
    colnames(umaps)=c('UMAP_1','UMAP_2')}
  if(dims==3){colnames(umaps)=c('UMAP_1','UMAP_2','UMAP_3')}
  if(dims==2){
    obj@reductions$UMAP<-as(umaps,'Matrix')
    return(obj)}else{return(umaps)}
}


### ****** 3. a function to run motif enrichment using ChromVAR standard pipeline ------------
### the input is a a binary (cell-by-peak) matrix and corresponding peaks' bed information !!!!!
RunChromVAR <- function(
  obj, ###a binary (cell-by-peak) matrix, the row of Bmat must = the number of peaks
  peak.obj = NULL, ### a GrangeList object of peak files or the directory path of ATAC peak files
  peakFormat = c('binary_matrix'), ## the format of peak object: GRangeList or peak bed file 'binary_matrix',,'peak_bed'
  Org=c('hg19'),#'mm10','hg38',
  min.count=1, ### the threshold of a peaks found at at least 10 cells
  species = c("Homo sapiens") ### default is Homo sapiens ,"Mus musculus"
){
  Bmat<-obj@bmat$filter
  options(stringsAsFactors = F)
  cat(">> checking depedent packages ... \t\t\t", format(Sys.time(), 
                                                         "%Y-%m-%d %X"), "\n")
  
  if (!(requireNamespace("SummarizedExperiment", quietly = TRUE))) {
    message("Please install package 'SummarizedExperiment'")
    BiocManager::install('SummarizedExperiment')
  }
  
  if (!(requireNamespace("chromVAR", quietly = TRUE))) {
    message("Please install package 'chromVAR'");
    BiocManager::install('chromVAR')
  }
  
  if (!(requireNamespace("GenomicFeatures", quietly = TRUE))) {
    cat("Please install package 'GenomicFeatures'");
    BiocManager::install('GenomicFeatures')
  }
  
  
  if (!(requireNamespace("GenomicRanges", quietly = TRUE))) {
    cat("Please install package 'GenomicRanges'");
    BiocManager::install('GenomicRanges')
  }
  
  if (!(requireNamespace("motifmatchr", quietly = TRUE))) {
    message("Please install package 'motifmatchr'")
    BiocManager::install('motifmatchr')
  }
  
  if (!(requireNamespace("JASPAR2016", quietly = TRUE))) {
    message("Please install package 'JASPAR2016'")
    BiocManager::install('JASPAR2016')
  }
  library(SummarizedExperiment)
  library(chromVAR)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(motifmatchr)
  library(JASPAR2016)
  cat(">>> checking input matrix ... \t\t\t", format(Sys.time(), 
                                                     "%Y-%m-%d %X"), "\n")
  
  
  if((x=max(Bmat)) > 1L){
    stop("input matrix is not binarized, run 'makeBinary' first!")
  }
  
  
  if(!is(min.count, "numeric")){
    stop("min.count is not a numeric");
  }else{
    min.count = max(min.count, 0);
  }
  
  #idy = which(Matrix::rowSums(Bmat) >= min.count);
  #Bmat = Bmat[,idy,dropping=TRUE]
  #peaksID = rownames(Bmat)
  
  if(dim(Bmat)[2] == 0){
    stop("input matrix is empty after filering low-coverage features")
  }
  
  cat(">>> checking genome annotation ... \t\t\t", format(Sys.time(), 
                                                          "%Y-%m-%d %X"), "\n")
  
  if(Org == 'mm9'){
    if (!(requireNamespace("BSgenome.Mmusculus.UCSC.mm9", quietly = TRUE))) {
      message("Please install package 'BSgenome.Mmusculus.UCSC.mm9'");
      BiocManager::install('BSgenome.Mmusculus.UCSC.mm9')
    }
    txdb = BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
    OrgDb = 'org.Mm.eg.db'
  }
  if(Org == 'mm10'){
    if (!(requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))) {
      message("Please install package 'BSgenome.Mmusculus.UCSC.mm10'");
      BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')
    }
    txdb = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    OrgDb = 'org.Mm.eg.db'
  }
  if(Org == 'hg19'){
    if (!(requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))) {
      message("Please install package 'BSgenome.Hsapiens.UCSC.hg19'");
      BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
    }
    txdb = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    OrgDb = 'org.Hs.eg.db'
  }
  if(Org == 'hg38'){
    if (!(requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))) {
      message("Please install package 'BSgenome.Hsapiens.UCSC.hg38'");
      BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
    }
    txdb = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    OrgDb = 'org.Hs.eg.db'
  }
  
  if(!(is(txdb, "BSgenome"))){
    stop("genome is not a BSgenome object");
  }
  
  
  if(!(is(species, "character"))){
    stop("species is not a character")
  }
  
  
  
  cat(">>> checking peaks ... \t\t\t", format(Sys.time(), 
                                              "%Y-%m-%d %X"), "\n")
  
  if(peakFormat=='binary_matrix'){
    # bins = as.character(rownames(Bmat))
    # bins_bed = do.call(rbind, strsplit(x = bins, split = '[:-]'))
    # bins_bed = as.data.frame(bins_bed);colnames(bins_bed)[1:3] = c('chr','start','end')
    # bins_bed[,2] = as.numeric(bins_bed[,2])
    # bins_bed[,3] = as.numeric(bins_bed[,3])
    # peak.use = GenomicRanges::GRanges(bins_bed[,1],IRanges::IRanges(start = bins_bed[,2],
    #                                                                 end = bins_bed[,3], names = bins))
    # peak.use = sort(peak.use)
    Bmat <- Bmat  
    bins = as.character(rownames(Bmat))
    bins_bed = do.call(rbind, strsplit(x = bins, split = '[:-]'))
    bins_bed = as.data.frame(bins_bed)
    
    write.table(bins_bed, "./bins_bed.bed", quote=F, sep="\t",row.names=F, col.names=F)
    peakfile <- "./bins_bed.bed"
    peak.use <- chromVAR::getPeaks(peakfile, sort_peaks = TRUE)
    file.remove("./bins_bed.bed")
  }else if (peakFormat == 'peak_bed'){
    peak.use <- chromVAR::getPeaks(peak.obj, sort_peaks = TRUE)
  }else if (sum(!is.null(peak.obj))>0){
    peak.use <- sort(peak.obj)
  }
  
  
  if((x=length(peak.use)) == 1L){
    stop("peak is empty!")
  }
  
  cat(">>> creating chromVAR object ...\t\t\t", format(Sys.time(), 
                                                       "%Y-%m-%d %X"), "\n")
  
  colData = S4Vectors::DataFrame(CellType = colnames(Bmat),
                                 depth=Matrix::colSums(Bmat))
  
  rse <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = Bmat), 
    rowRanges = peak.use, 
    colData = colData)
  
  cat(">>> computing GC bias then filtering... \t\t\t", format(Sys.time(), 
                                                               "%Y-%m-%d %X"), "\n")
  set.seed(2020)
  rse <- chromVAR::addGCBias(rse, genome = txdb);
  row.data <- data.frame(rowData(rse))
  row.data[is.na(row.data)] <- 0
  rowData(rse) <- row.data
  counts_filtered <- chromVAR::filterPeaks(rse,min_fragments_per_peak=min.count)
  
  cat(">>> getting JASPAR motifs ...\t\t\t", format(Sys.time(), 
                                                    "%Y-%m-%d %X"), "\n")
  motifs <- chromVAR::getJasparMotifs(collection = "CORE", species = species);
  
  
  cat(">>> scanning motifs in the peaks ... \t\t\t", format(Sys.time(), 
                                                            "%Y-%m-%d %X"), "\n")
  motif_mm <- motifmatchr::matchMotifs(motifs, counts_filtered, genome = txdb);
  
  cat(">>> calculating motif variability between cells ... \t\t\t", format(Sys.time(), 
                                                                           "%Y-%m-%d %X"), "\n")
  dev <- chromVAR::computeDeviations(object = counts_filtered, annotations = motif_mm);
  dev_mat = t(SummarizedExperiment::assay(dev));
  
  
  cat(">>> Done ... \t\t\t", format(Sys.time(), 
                                    "%Y-%m-%d %X"), "\n")
  obj@mmat<-t(as(dev_mat,'Matrix'))
  return(obj)
}

### peaks were directly extracted from the rownames of binary matrix
## hg19_motif = RunChromVAR(peak.obj = NULL,peakFormat = 'binary_matrix', 
##                         Bmat = binary_matrix_filtered,
##                         Org='hg19', species="Homo sapiens")

### manually prepare a peak GrangeList
## hg19_motif = RunChromVAR(peak.obj =bins.gr,peakFormat = 'GRangeList', 
##                         Bmat = binary_matrix_filtered,
##                         Org='hg19', species="Homo sapiens")



### ****** 2. a function to visualize the ATAC signals (gene levels) in UMAP/tSNE embeddings  --------
PlotSelectGenesATAC = function(obj, gene2plot = c("Snap25", "Gad2", "Apoe",'BCL9'), 
                               reduction = 'TSNE',
                               ncol = NULL){
  plot_theme <- theme(plot.title = element_text(hjust = 0.5, size = 20),
                      #legend.position = 'right',
                      legend.title =element_text(size=15),
                      legend.text = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      axis.title.x = element_text(size=15),
                      axis.title.y = element_text(size=15),
                      axis.text.y  = element_text(size=15),
                      panel.border = element_blank(),
                      axis.line.x = element_line(size=0.25, color="black"),
                      axis.line.y = element_line(size=0.25, color="black"),
                      panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
                      panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                      panel.background = element_rect(fill='white'),
                      legend.key=element_blank(),
                      strip.text.x = element_text(size=15),
                      strip.text.y = element_text(size=15),
                      strip.background = element_rect(colour = 'white', fill = 'white'))
  
  if (!(requireNamespace("ggplot2", quietly = TRUE))) {
    cat("Please install package 'ggplot2'");
    install.packages('ggplot2')
  }
  
  if (!(requireNamespace("ggpubr", quietly = TRUE))) {
    cat("Please install package 'ggpubr'");
    install.packages('ggpubr')
  }
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  Gmat<-obj@gmat
  cellsReduction <- data.frame(obj@barcode)
  
  # if(sum(!is.data.frame(Gmat))>0){
  #   Gmat = as.data.frame(Gmat)
  # }
  lapply(gene2plot, function(i){if(!(i %in% rownames(Gmat) ))
  {message(paste0(i, ' doesn not exist' ))}} )
  if (sum(gene2plot%in%rownames(Gmat))==0){stop('no gene is found')}
  mat2plot =  Gmat[rownames(Gmat)%in%gene2plot,]
  if(mode(mat2plot)=='numeric'){ 
    mat2plot=as.data.frame(t(mat2plot))

    rownames(mat2plot)<-rownames(Gmat)[rownames(Gmat)%in%gene2plot]}
  gene2plot = rownames(mat2plot)
  
  if(length(gene2plot)>9){
    stop('The maximun of genes to plot was limmited to 9!')
  }
  p.list = lapply(1:nrow(mat2plot),
                  function(i){
                    gene = rownames(mat2plot)[i]
                    cellsReduction[,'level'] = as.numeric(mat2plot[gene,])
                    od <- order(cellsReduction[,'level'])
                    cellsReduction <- cellsReduction[od,]
                    if(reduction == 'UMAP'){
                      UMAP_1<-obj@reductions$UMAP[od,'UMAP_1']
                      UMAP_2<-obj@reductions$UMAP[od,'UMAP_2']
                      umap = ggplot(cellsReduction[,],
                                    aes(x=UMAP_1,y=UMAP_2,col=level))+ geom_point(size=1) +
                        scale_color_gradient(low = 'lightgrey', high = 'red') +
                        labs(title = paste("Counts, ",gene, sep=""), x = 'UMAP_1', y = 'UMAP_2')+  
                        plot_theme + theme(legend.position = 'right', legend.title = element_blank()) 
                      umap
                    }else{
                      tSNE_1<-obj@reductions$TSNE@matrix[od,'tSNE_1']
                      tSNE_2<-obj@reductions$TSNE@matrix[od,'tSNE_2']
                      tsne = ggplot(cellsReduction[,],
                                    aes(x=tSNE_1,y=tSNE_2,col=level))+ geom_point(size=1) +
                        # scale_color_gradient(low = 'yellow', high = 'blue') +
                        scale_color_gradient(low = 'lightgrey', high = 'red') +
                        labs(title = paste("Counts, ",gene, sep=""), x = 'tSNE_1', y = 'tSNE_2')+
                        plot_theme + theme(legend.position = 'right', legend.title = element_blank())
                      tsne
                    }
                  })
  
  
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = gene2plot) == 1) {
      ncol <- 1
    }
    if (length(x = gene2plot) > 4) {
      ncol <- 3
    }
  }
  
  if(length(gene2plot)==1){
    p=p.list[[1]]
    print(p.list[[1]])
  }
  if(length(gene2plot)==2){
    p=plot_grid(p.list[[1]], p.list[[2]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]], ncol=ncol))
  }
  if(length(gene2plot)==3){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]], ncol=ncol))
  }
  if(length(gene2plot)==4){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]], ncol=ncol))
  }
  if(length(gene2plot)==5){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]], p.list[[5]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]], p.list[[5]], ncol=ncol))
  }
  if(length(gene2plot)==6){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]],
                p.list[[5]], p.list[[6]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]],
                    p.list[[5]], p.list[[6]], ncol=ncol))
  }
  if(length(gene2plot)==7){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]],
                p.list[[5]], p.list[[6]],
                p.list[[7]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]],
                    p.list[[5]], p.list[[6]],
                    p.list[[7]], ncol=ncol))
  }
  if(length(gene2plot)==8){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]],
                p.list[[5]], p.list[[6]],
                p.list[[7]], p.list[[8]],ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]],
                    p.list[[5]], p.list[[6]],
                    p.list[[7]], p.list[[8]],ncol=ncol))
  }
  if(length(gene2plot)==9){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]],
                p.list[[5]], p.list[[6]],
                p.list[[7]], p.list[[8]],
                p.list[[9]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]],
                    p.list[[5]], p.list[[6]],
                    p.list[[7]], p.list[[8]],
                    p.list[[9]], ncol=ncol))
  }
  return(p)
}
mcesApply <- function(X, MARGIN, FUN, cores=1, ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  multiassign(names(pData(X)), pData(X), envir=e1)
  environment(FUN) <- e1
  cl <- makeCluster(cores)
  clusterEvalQ(cl, {require(VGAM);})
  if (MARGIN == 1){
    res <- parRapply(cl, exprs(X), FUN, ...)
  }else{
    res <- parCapply(cl, exprs(X), FUN, ...)
  }
  
  stopCluster(cl)
  res
}


diff_test_helperBeta <- function(x, 
                                 fullModelFormulaStr, 
                                 reducedModelFormulaStr, 
                                 relative_expr,
                                 weights,
                                 disp_func=NULL,
                                 verbose=FALSE
){ 
  
  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, sep="")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, sep="")
  
  x_orig <- x
  disp_guess <- 0
  
  # else if (expressionFamily@vfamily %in% c("binomialff")){
  f_expression <- x/Size_Factor
  #f_expression[f_expression > 1] <- 1
  # }
  
  
  test_res <- tryCatch({
    # if (expressionFamily@vfamily %in% c("binomialff")){
    #   if (verbose){
    expressionFamily="binomialff"
    full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
    reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)                         
    #   }else{
    #     full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
    #     reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))                    
    #   }
    # }
    
    #print(full_model_fit)
    #print(coef(reduced_model_fit))
    compareModelsBeta <- function(full_models, reduced_models){
      stopifnot(length(full_models) == length(reduced_models))
      test_res <- mapply(function(x,y) { 
        if (is.null(x) == FALSE && is.null(y) == FALSE) {
          lrt <- VGAM::lrtest(x,y) 
          pval=lrt@Body["Pr(>Chisq)"][2,]
          family = x@family@vfamily
          if (length(family) > 1)
            family = family[1]
          beta = x@coefficients[2]
          data.frame(status = "OK", family=family, pval=pval,beta=beta)
        } else { data.frame(status = "FAIL", family=NA, pval=1.0,beta=0) } 
      } , full_models, reduced_models, SIMPLIFY=FALSE, USE.NAMES=TRUE)
      
      test_res <- do.call(rbind.data.frame, test_res)
      test_res$qval <- p.adjust(test_res$pval, method="BH")
      test_res
    }
    
    compareModelsBeta(list(full_model_fit), list(reduced_model_fit))
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { 
    if(verbose)
      print (e);
    data.frame(status = "FAIL", pval=1.0, qval=1.0, beta=0)
    #data.frame(status = "FAIL", pval=1.0) 
  }
  )
  test_res
}

differentialGeneTest <- function(cds, 
                                 fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
                                 reducedModelFormulaStr="~1", 
                                 relative_expr=TRUE,
                                 cores=1, 
                                 verbose=FALSE
){
  status <- NA
  
  if (cores > 1){
    diff_test_res<-mcesApply(cds, 1, diff_test_helperBeta, 
                             c("BiocGenerics", "VGAM", "Matrix"), 
                             cores=cores, 
                             fullModelFormulaStr=fullModelFormulaStr,
                             reducedModelFormulaStr=reducedModelFormulaStr,
                             relative_expr=relative_expr,
                             disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                             verbose=verbose
                             
    )
    diff_test_res
  }else{
    diff_test_res<-smartEsApply(cds,1,diff_test_helperBeta, 
                                convert_to_dense=TRUE,
                                fullModelFormulaStr=fullModelFormulaStr,
                                reducedModelFormulaStr=reducedModelFormulaStr, 
                                relative_expr=relative_expr,
                                disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                                verbose=verbose
                                #          ,
                                # backup_method = backup_method, 
                                # use_epislon = use_epislon,
                                # stepsize = stepsize
                                
    )
    diff_test_res
  }
  
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)
  
  diff_test_res$qval <- 1
  diff_test_res$qval[which(diff_test_res$status == 'OK')] <- p.adjust(subset(diff_test_res, status == 'OK')[, 'pval'], method="BH")
  
  diff_test_res <- merge(diff_test_res, fData(cds), by="row.names")
  row.names(diff_test_res) <- diff_test_res[, 1] #remove the first column and set the row names to the first column
  diff_test_res[, 1] <- NULL 
  
  diff_test_res
}
# test_clust=1
# data<-art@bmat$filter[1:5,]
# cluster_result<-art@metaData[1,5]



FindDAR <- function(obj, test_clust,out_dir=NULL){
  if(is.null(out_dir)){out_dir=getwd()}
  cluster <- test_clust
  Cluster <- paste("Cluster",cluster,sep="")
  library(Matrix)
  da_mat_final = as.matrix(obj@bmat$filter)
  dclust <- obj@metaData$cluster
  clusterfactor <- factor((dclust==as.character(cluster))+0,levels=c("0","1"))
  pda = data.frame(as.character(colnames(da_mat_final)),
                   clusterfactor)
  names(pda) = c("CellID","CellCluster")
  pda = new("AnnotatedDataFrame", data = pda)
  fda = as.data.frame(as.character(rownames(da_mat_final)))
  names(fda) = "Peak"
  fda = new("AnnotatedDataFrame", data = fda)
  rownames(da_mat_final) = NULL
  colnames(da_mat_final) = NULL
  library(SummarizedExperiment)
  library(monocle)
  submat_cds = suppressWarnings( newCellDataSet(da_mat_final,
                                                featureData = fda,
                                                phenoData = pda,
                                                lowerDetectionLimit=1))
  
  pData(submat_cds)$Size_Factor = 1
  differtest =  suppressWarnings( differentialGeneTest(submat_cds, fullModelFormulaStr = "~CellCluster ",cores=15))
  write.table(differtest[order(differtest$qval),], paste(out_dir,'/',Cluster,"_peak.txt",sep=""),quote=F,sep="\t")
  diff <- differtest[order(differtest$qval),]
  diff <- diff[diff[,"qval"]<0.001,]
  write.table(diff, paste(Cluster,"_specific_peak.txt",sep=""),quote=F,sep="\t")
  dt <- as.data.frame(diff[,"Peak"])
  parts <-function(x) {
    m <- regexec("([a-z,A-Z | 0-9]+)(-|,)*([0-9]*)(-|,)*([0-9]*)", x)
    parts <- unlist(lapply(regmatches(x, m), '[',  c(2L, 4L, 6L)))
    parts
  }
  
  rr<-do.call(rbind,lapply(as.character(dt[,1]),parts))
  write.table(rr, paste("./",Cluster,"_specific_peak.bed",sep=""),quote=F,sep="\t",row.names=F, col.names=F)
  
  return(differtest)
}






plotTrajectory <- function(obj,anno=NULL,color=NULL) {
  if (is.null(anno)) {
    anno <- obj@metaData$cluster
  } else {
    anno <- eval(parse(text =paste0('obj@metaData$',anno) ))
  }
  library(RColorBrewer)
  if (is.null(color)) {
    tsnecols = c(brewer.pal(12, "Paired"), brewer.pal(7, "Accent"), brewer.pal(9, "Pastel1"))
    tsnecols = unique(tsnecols)
  } else {
    tsnecols = color
  }
  library(scatterplot3d)
  
  points <- obj@trajectory
  mydata<-cbind(points[,1],points[,2],points[,3])
  colnames(mydata)<-c("1m","2m","3m")
  mydata<-as.data.frame(mydata)
  with(mydata, {
    scatterplot3d(points[,3],   # x axis
                  points[,2],     # y axis
                  points[,1], 
                  color=tsnecols[as.factor(anno)],
                  pch=19, cex.symbols=0.4,xlab="dim 1",ylab="dim 2",zlab="dim 3",cex.axis=0.75)
  })
  legend("topright",pch=c(rep(20,length(names(table(anno))))),col=tsnecols[1:length(table(as.factor(anno)))],
         legend=names(table(anno)))
  
  
}


RunTrajectory <- function(obj, anno=NULL,sigma=NULL, lambda=NULL, nSV=NULL, ndim=NULL,data=NULL,gamma=NULL){
  
  if (is.null(sigma)) {
    sigma <- 0.001
  } else {
    sigma <- sigma
  }
  if (is.null(lambda)) {
    lambda <- NULL
  } else {
    lambda <- lambda
  }
  if (is.null(anno)) {
    anno <- as.data.frame(obj@metaData$cluster)
    rownames(anno)<-rownames(obj@metaData)
  } else {
    anno <- as.data.frame(eval(parse(text =paste0('obj@metaData$',anno) )))
    rownames(anno)<-rownames(obj@metaData)
  }
  if (is.null(nSV)) {
    nSV <- 10
  } else {
    nSV <- nSV
  }
  if (is.null(ndim)) {
    ndim <- 5
  } else {
    ndim <- ndim
  }
  if (is.null(gamma)) {
    gamma <- 10
  } else {
    gamma <- gamma
  }
  
  decomp_data  <- obj@reductions$SVD@x 
  pca.var <- obj@reductions$SVD@sdev
  
  cut <- pca.var[1]/pca.var[2]
  if (cut > 2){
    svd_tsne = decomp_data[,2:nSV]
  }else{
    svd_tsne = decomp_data[,1:nSV]
  }
  print(dim(svd_tsne)[2])
  
  library(Rtsne)
  set.seed(10)
  data <- svd_tsne
  tsne_out = Rtsne(data, pca=F, dims =5 , perplexity = 30,
                   theta = 0.5, check_duplicates = TRUE, max_iter = 500)
  data <- tsne_out$Y
  
  cellnames <- obj@barcode
  library(monocle)
  expr_matrix = t(data)
  colnames(expr_matrix) <- obj@barcode
  rownames(expr_matrix) <- c(paste("dim",c(1:dim(expr_matrix)[1]),sep="_"))
  
  cells <- as.data.frame(anno)
  
  genes <- as.data.frame(rownames(expr_matrix))
  colnames(genes) <- c("gene_short_name")
  rownames(genes) <- rownames(expr_matrix)
  
  pd <- new("AnnotatedDataFrame", data=cells)
  fd <- new("AnnotatedDataFrame", data=genes)
  HSMM <- newCellDataSet(as.matrix(expr_matrix), 
                         phenoData=pd, featureData=fd, 
                         expressionFamily=negbinomial.size())
  HSMM <- suppressWarnings(estimateSizeFactors(HSMM))
  
  HSMM <- suppressWarnings(reduceDimension(HSMM, max_components=3, norm_method='none', 
                                           scaling=T,reduction_method="DDRTree", 
                                           sigma = sigma,lambda=lambda, param.gamma = gamma, verbose=F))
  
  points <- t(HSMM@reducedDimS)
  obj@trajectory<- as(points,'Matrix')
  return(obj)
}




MapBin2Gene = function(obj, ### the cell-by-bin matrix
                       binFormat = 'binary_matrix', ### the format of cell-by-bin matrix
                       bin_file = NULL,
                       Org = 'hg19', ### mm10
                       OrgDb = 'org.Hs.eg.db', ### 
                       TxDb = NULL, ### if Org = manual, you should input the TxDb defined by yourself 
                       convert_mat = TRUE, ### whether convert bin-by-cell matrix to cell-by
                       TSS_window = 3000 ### the window size around TSS to define the promoter 
){
  
  options(stringsAsFactors = F)
  
  cat(">> checking depedent packages ... \t\t\t", format(Sys.time(), 
                                                         "%Y-%m-%d %X"), "\n")
  Bmat<-obj@bmat$filter
  if (!(requireNamespace("BiocManager", quietly = TRUE))) {
    cat("Please install package 'BiocManager'");
    install.packages('BiocManager')
  }
  
  if (!(requireNamespace("data.table", quietly = TRUE))) {
    cat("Please install package 'data.table'");
    install.packages('data.table')
  }
  
  if (!(requireNamespace("GenomicFeatures", quietly = TRUE))) {
    cat("Please install package 'GenomicFeatures'");
    BiocManager::install('GenomicFeatures')
  }
  
  if (!(requireNamespace("ChIPseeker", quietly = TRUE))) {
    cat("Please install package 'ChIPseeker'");
    BiocManager::install('ChIPseeker')
  }
  
  
  
  ### load the ATAC IDs of peaks or bins -------  
  cat(">> loading ATAC features ...\t\t\t", format(Sys.time(), 
                                                   "%Y-%m-%d %X"), "\n")
  
  if(sum(!is.matrix(Bmat))>0){
    Bmat = Matrix::as.matrix(Bmat)
    ### matrix is faster than data.frame
  }
  
  if(binFormat == 'binary_matrix'){
    bins = as.character(rownames(Bmat))
    bins_bed = do.call(rbind, strsplit(x = bins, split = '[:-]'))
    bins_bed = as.data.frame(bins_bed);
    colnames(bins_bed)[1:3] = c('chr','start','end')
    bins_bed[,2] = as.numeric(bins_bed[,2])
    bins_bed[,3] = as.numeric(bins_bed[,3])
    bins.gr = GenomicRanges::GRanges(bins_bed[,1], IRanges(start = bins_bed[,2], end = bins_bed[,3], names = bins))
  }
  if(binFormat == 'bin_file'){
    bins_bed = read.table(bin_file, sep='\t', fill =T, stringsAsFactors = F, header = F)
    colnames(bins_bed)[1:3] = c('chr','start','end')
    bins_bed = transform(bins_bed, name = paste0(chr,':',start,'-',end))
    bins.gr = GenomicRanges::GRanges(bins_bed$chr, IRanges(start = bins_bed$start, end = bins_bed$end, names = bins))
  }
  
  
  ### load the genome reference  -------- 
  cat(">> adding gene annotation...\t\t\t", format(Sys.time(), 
                                                   "%Y-%m-%d %X"), "\n")
  
  if(Org == 'mm9'){
    if (!(requireNamespace("TxDb.Mmusculus.UCSC.mm9.knownGene", quietly = TRUE))) {
      message("Please install package 'TxDb.Mmusculus.UCSC.mm9.knownGene'");
      BiocManager::install('TxDb.Mmusculus.UCSC.mm9.knownGene')
    }
    txdb = TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
    OrgDb = 'org.Mm.eg.db'
  }
  if(Org == 'mm10'){
    if (!(requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE))) {
      message("Please install package 'TxDb.Mmusculus.UCSC.mm10.knownGene'");
      BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene')
    }
    txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    OrgDb = 'org.Mm.eg.db'
  }
  if(Org == 'hg19'){
    if (!(requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE))) {
      message("Please install package 'TxDb.Hsapiens.UCSC.hg19.knownGene'");
      BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')
    }
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    OrgDb = 'org.Hs.eg.db'
  }
  if(Org == 'hg38'){
    if (!(requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))) {
      message("Please install package 'TxDb.Hsapiens.UCSC.hg38.knownGene'");
      BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')
    }
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    OrgDb = 'org.Hs.eg.db'
  }
  if(Org == 'manual'){
    txdb = TxDb
  }
  
  ### map to genome reference -----
  cat(">> mapping bins to genome...\t\t\t", format(Sys.time(), 
                                                   "%Y-%m-%d %X"), "\n")
  
  binAnno <- ChIPseeker::annotatePeak(bins.gr, tssRegion=c(-TSS_window, TSS_window),
                                      TxDb=txdb, annoDb=OrgDb, level = 'gene',
                                      genomicAnnotationPriority = c("Promoter","Exon", "Intron",
                                                                    "Downstream", "Intergenic")) 
  
  bin2gene = as.data.frame(binAnno@anno)
  library(data.table)
  data.table::setDT(bin2gene,keep.rownames = 'ID')
  head(bin2gene)
  
  bin2gene[,Location:=annotation]
  bin2gene[grep('Exon',annotation)]$Location = 'Exon'
  bin2gene[grep('Intron',annotation)]$Location = 'Intron'
  bin2gene[grep('Promoter',annotation)]$Location = 'Promoter'
  bin2gene[grep('Down',annotation)]$Location = 'Downstream'
  
  bin2use = bin2gene[grep("Promoter",Location)]
  gene2use = unique(bin2use$SYMBOL)
  
  ### select the ATAC bins in gene's promoter
  
  ### convert cell-by-bin matrix to cell-by-gene matrix 
  mat2plot = do.call(rbind, lapply(gene2use, function(gene){
    bins2use = bin2use[SYMBOL==gene]$ID
    if(length(bins2use)==1){
      sum_bins = Bmat[bins2use,]
    }else{
      sum_bins = colSums(Bmat[bins2use,])
    }
    sum_bins
  }))
  rownames(mat2plot) = gene2use
  cat(">>", nrow(mat2plot), "genes with promoter ATAC bins retained after filtering ...\t\t\t", format(Sys.time(),  "%Y-%m-%d %X"), "\n")  
  # mat2plot<- t(t(mat2plot)/art@metaData$nCounts)*1000000
  obj@gmat<-as(as.matrix(mat2plot), "dgCMatrix")
  return(obj)
}

PlotSelectTF= function(obj,TF2plot = c("GSC2", "EVX1566", "GSX2"), 
                       reduction = 'TSNE',
                       ncol = NULL){
  plot_theme <- theme(plot.title = element_text(hjust = 0.5, size = 20),
                      #legend.position = 'right',
                      legend.title =element_text(size=15),
                      legend.text = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      axis.title.x = element_text(size=15),
                      axis.title.y = element_text(size=15),
                      axis.text.y  = element_text(size=15),
                      panel.border = element_blank(),
                      axis.line.x = element_line(size=0.25, color="black"),
                      axis.line.y = element_line(size=0.25, color="black"),
                      panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
                      panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                      panel.background = element_rect(fill='white'),
                      legend.key=element_blank(),
                      strip.text.x = element_text(size=15),
                      strip.text.y = element_text(size=15),
                      strip.background = element_rect(colour = 'white', fill = 'white'))
  
  if (!(requireNamespace("ggplot2", quietly = TRUE))) {
    cat("Please install package 'ggplot2'");
    install.packages('ggplot2')
  }
  
  if (!(requireNamespace("ggpubr", quietly = TRUE))) {
    cat("Please install package 'ggpubr'");
    install.packages('ggpubr')
  }
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  Gmat<-obj@mmat
  cellsReduction <- data.frame(obj@barcode)
  Gmat <- Gmat[,obj@barcode]
  # if(sum(!is.data.frame(Gmat))>0){
  #   Gmat = as.data.frame(Gmat)
  # }
  
  lapply(TF2plot, function(i){
    if(length(grep(i,rownames(Gmat)))==0)
    {message(paste0(i, ' doesn not exist' ))}
  } )
  if (length(grep(paste(TF2plot,collapse = '|'),rownames(Gmat)))==0)
    {stop('no TF is found')}
  mat2plot =  Gmat[grep(paste(TF2plot,collapse = '|'),rownames(Gmat)),]
  
  if(mode(mat2plot)=='numeric'){ 
    mat2plot=as.data.frame(t(mat2plot))
    rownames(mat2plot)<-rownames(Gmat)[grep(paste(TF2plot,collapse = '|'),rownames(Gmat))]}
  TF2plot = rownames(mat2plot)
  
  if(length(TF2plot)>9){
    stop('The maximun of genes to plot was limmited to 9!')
  }
  
  p.list = lapply(1:nrow(mat2plot),
                  function(i){
                    gene = rownames(mat2plot)[i]
                    cellsReduction[,'level'] = as.numeric(mat2plot[gene,])
                    # od <- order(cellsReduction[,'level'])
                    # cellsReduction <- cellsReduction[od, ]
                    if(reduction == 'UMAP'){
                      UMAP_1<-obj@reductions$UMAP[,'UMAP_1']
                      UMAP_2<-obj@reductions$UMAP[,'UMAP_2']
                      umap = ggplot(cellsReduction,
                                    aes(x=UMAP_1,y=UMAP_2,col=level))+ geom_point(size=1) +
                        scale_color_gradient2(low = 'blue', high = 'red') +
                        labs(title = gene, x = 'UMAP_1', y = 'UMAP_2')+  
                        plot_theme +theme(legend.position = 'right', legend.title = element_blank())
                      umap
                    }else{
                      tSNE_1<-obj@reductions$TSNE@matrix[,'tSNE_1']
                      tSNE_2<-obj@reductions$TSNE@matrix[,'tSNE_2']
                      tsne = ggplot(cellsReduction,
                                    aes(x=tSNE_1,y=tSNE_2,col=level))+ geom_point(size=1) +
                        scale_color_gradient2(low = 'blue', high = 'red') +
                        labs(title = gene, x = 'tSNE_1', y = 'tSNE_2')+
                        plot_theme +theme(legend.position = 'right', legend.title = element_blank())
                      tsne
                    }
                  })
  
  
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = TF2plot) == 1) {
      ncol <- 1
    }
    if (length(x = TF2plot) > 4) {
      ncol <- 3
    }
  }
  
  if(length(TF2plot)==1){
    p=p.list[[1]]
    print(p.list[[1]])
  }
  if(length(TF2plot)==2){
    p=plot_grid(p.list[[1]], p.list[[2]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]], ncol=ncol))
  }
  if(length(TF2plot)==3){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]], ncol=ncol))
  }
  if(length(TF2plot)==4){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]], ncol=ncol))
  }
  if(length(TF2plot)==5){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]], p.list[[5]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]], p.list[[5]], ncol=ncol))
  }
  if(length(TF2plot)==6){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]],
                p.list[[5]], p.list[[6]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]],
                    p.list[[5]], p.list[[6]], ncol=ncol))
  }
  if(length(TF2plot)==7){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]],
                p.list[[5]], p.list[[6]],
                p.list[[7]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]],
                    p.list[[5]], p.list[[6]],
                    p.list[[7]], ncol=ncol))
  }
  if(length(TF2plot)==8){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]],
                p.list[[5]], p.list[[6]],
                p.list[[7]], p.list[[8]],ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]],
                    p.list[[5]], p.list[[6]],
                    p.list[[7]], p.list[[8]],ncol=ncol))
  }
  if(length(TF2plot)==9){
    p=plot_grid(p.list[[1]], p.list[[2]],
                p.list[[3]],p.list[[4]],
                p.list[[5]], p.list[[6]],
                p.list[[7]], p.list[[8]],
                p.list[[9]], ncol=ncol)
    print(plot_grid(p.list[[1]], p.list[[2]],
                    p.list[[3]],p.list[[4]],
                    p.list[[5]], p.list[[6]],
                    p.list[[7]], p.list[[8]],
                    p.list[[9]], ncol=ncol))
  }
  return(p)
}


# Read_10X <- function(
#   data.dir,
#   cell.column = 1,
#   unique.features = TRUE,
#   strip.suffix = FALSE) {
#   for (i in seq_along(along.with = data.dir)) {
#     run <- data.dir[i]
#     if (!dir.exists(paths = run)) {
#       stop("Directory provided does not exist")
#     }
#     barcode.loc <- file.path(run, 'barcodes.tsv')
#     peaks.loc <- file.path(run, 'peaks.bed')
#     matrix.loc <- file.path(run, 'matrix.mtx')
#     # Flag to indicate if this data is from CellRanger >= 3.0
#     
#     if (!file.exists(barcode.loc)) {
#       stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
#     }
#     if (!file.exists(peaks.loc) ) {
#       stop("Peaks name or features file missing. Expecting ", basename(path = features.loc))
#     }
#     if (!file.exists(matrix.loc)) {
#       stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
#     }
#     library(Matrix)
#     data <- as(readMM(file = matrix.loc),'dgCMatrix')
#     cell.barcodes <- read.table(file = barcode.loc, header = FALSE, row.names = NULL)
#     
#     
#     peaks <- read.delim(
#       file = peaks.loc,
#       header = FALSE, stringsAsFactors = FALSE)
#     
#     colnames(data)<-cell.barcodes$V1
#     rownames(data)<-paste0(peaks$V1,'-',peaks$V2,'-',peaks$V3)
#     art<-CreatescART(data)
#     return(art)
#   }
# }
# 
# 
# Read_10X_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
#   if (!requireNamespace('hdf5r', quietly = TRUE)) {
#     stop("Please install hdf5r to read HDF5 files")
#   }
#   if (!file.exists(filename)) {
#     stop("File not found")
#   }
#   infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
#   genomes <- names(x = infile)
#   output <- list()
#   if (hdf5r::existsGroup(infile, 'matrix')) {
#     # cellranger version 3
#     if (use.names) {
#       feature_slot <- 'features/name'
#     } else {
#       feature_slot <- 'features/id'
#     }
#   } else {
#     if (use.names) {
#       feature_slot <- 'gene_names'
#     } else {
#       feature_slot <- 'genes'
#     }
#   }
#   for (genome in genomes) {
#     counts <- infile[[paste0(genome, '/data')]]
#     indices <- infile[[paste0(genome, '/indices')]]
#     indptr <- infile[[paste0(genome, '/indptr')]]
#     shp <- infile[[paste0(genome, '/shape')]]
#     features <- infile[[paste0(genome, '/', feature_slot)]][]
#     barcodes <- infile[[paste0(genome, '/barcodes')]]
#     sparse.mat <- sparseMatrix(
#       i = indices[] + 1,
#       p = indptr[],
#       x = as.numeric(x = counts[]),
#       dims = shp[],
#       giveCsparse = FALSE
#     )
#     if (unique.features) {
#       features <- make.unique(names = features)
#     }
#     rownames(x = sparse.mat) <- features
#     colnames(x = sparse.mat) <- barcodes[]
#     sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
#     # Split v3 multimodal
#     if (infile$exists(name = paste0(genome, '/features'))) {
#       types <- infile[[paste0(genome, '/features/feature_type')]][]
#       types.unique <- unique(x = types)
#       if (length(x = types.unique) > 1) {
#         message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
#         sparse.mat <- sapply(
#           X = types.unique,
#           FUN = function(x) {
#             return(sparse.mat[which(x = types == x), ])
#           },
#           simplify = FALSE,
#           USE.NAMES = TRUE
#         )
#       }
#     }
#     output[[genome]] <- sparse.mat
#   }
#   infile$close_all()
#   output$matrix@Dimnames[[1]]=gsub(':','-', output$matrix@Dimnames[[1]])
#   obj<-CreatescART(output$matrix)
#   return(obj)
# }

Read_counts <- function(
  data.dir,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE) {
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, 'barcodes.tsv')
    peaks.loc <- file.path(run, 'bins.bed')
    matrix.loc <- file.path(run, 'matrix.mtx')
    # Flag to indicate if this data is from CellRanger >= 3.0
    
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!file.exists(peaks.loc) ) {
      stop("peaks name or features file missing. Expecting ", basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
    }
    library(Matrix)
    data <- as(readMM(file = matrix.loc),'dgCMatrix')
    cell.barcodes <- read.table(file = barcode.loc, header = FALSE, row.names = NULL)
    
    
    peaks <- read.delim(
      file = peaks.loc,
      header = FALSE, stringsAsFactors = FALSE)
    
    colnames(data)<-cell.barcodes$V1
    rownames(data)<-paste0(peaks$V1,'-',peaks$V2,'-',peaks$V3)
    art<-CreatescART(data)
    return(art)
  }
}

Snap2art<-function(snap){
  counts=t(snap@bmat)
  barcode=snap@barcode
  meta=snap@metaData
  pmat=snap@pmat
  bins=snap@feature
  art=CreatescART(data = counts,barcode = barcode,bins = bins)
  art@pmat=snap@pmat
  return(art)
}

Read_snap<-function(file,barcode,bin.size=NULL,sample="atac"){
  if (is.null(bin.size)) {
    bin.size <- 5000
  } else {
    bin.size  <- bin.size 
  }
  
  library(SnapATAC)
  sample <- sample
  x.sp <- createSnap(file = file,sample = sample)
barcodes = read.csv(
    barcode,
    head=TRUE
  );
barcodes = barcodes[2:nrow(barcodes),];
promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
UMI = log(barcodes$passed_filters+1, 10);
data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio);
barcodes$promoter_ratio = promoter_ratio;
library(viridisLite);
library(ggplot2);

barcodes.sel = barcodes[which(UMI >= 3 & UMI <= 5 & promoter_ratio >= 0.15 & promoter_ratio <= 0.6),];
rownames(barcodes.sel) = barcodes.sel$barcode;
x.sp = x.sp[which(x.sp@barcode %in% barcodes.sel$barcode),];
x.sp@metaData = barcodes.sel[x.sp@barcode,];
x.sp
bin.size <- 5000
x.sp = addBmatToSnap(x.sp, bin.size=bin.size, num.cores=1)
x.sp = makeBinary(x.sp, mat="bmat")
data <- t(x.sp@bmat)
length(x.sp@barcode)
length(x.sp@feature)
dim(data)
colnames(data) <- x.sp@barcode
rownames(data) <- x.sp@feature$name
head(colnames(data))

  obj<-new('scART')
  obj@metaData <- data.frame(x.sp@barcode)
  rownames(obj@metaData) <- x.sp@barcode
 
  obj@bmat=list(NULL,NULL,NULL,NULL,NULL)
  names(obj@bmat)=c('raw','binary','imputation','filter','TF_IDF')
  obj@bmat$raw=data
  obj@bmat$bmat=data
  obj@barcode<-x.sp@barcode
  obj@feature<-x.sp@feature
  nCounts<-sparse_Sums(data, rowSums = F)
  obj@metaData$nCounts<-nCounts
  return(obj)
}






#' Convert a snap object to a seurat object
#'
#' @param obj A snap object.
#' @param eigs.dims A vector of the dimensions to use.
#' @param norm A logical variable indicates whether to normalize the cell by gene matrix.
#' @param scale A logical variable indicagtes whether to scale the cell by gene matrix.
#' @return Return a seurat object that contains single cell ATAC-seq data
#'
#' @export

art2seurat <- function(
  obj, 
  eigs.dims=1:20,
  input.mat ="bmat",
  norm=TRUE,
  scale=TRUE){
  cat("Epoch: checking input parameters ... \n", file = stderr())
  if(missing(obj)){
    stop("obj is missing")
  }else{
    if(!is(obj,'scART')){
      stop("obj is not a snap object");
    }  
    if((x=length(obj@barcode))==0L){
      stop("obj@barcode is empty");   
    }
  }
  # check if Seurat is installed
  if (requireNamespace("Seurat", quietly = TRUE)) {
    require(Seurat)
  } else {
    stop("Please install Seurat V3 - learn more at https://github.com/satijalab/seurat");
  }
  
  if((x=nrow(obj@gmat)) == 0L){
    stop("gmat in obj is empty!");
  }
  gmat.use = (as(obj@gmat,'dgCMatrix'))
  # check if mat is binary;
  
  
  if(input.mat == "bmat"){
    data.use = obj@bmat$filter
    peak.use = as.data.frame(obj@feature);
  }else{
    data.use = obj@gmat;
    peak.use = as.data.frame(rownames(data.use));
  }
  
  if((x=nrow(data.use)) == 0L){
    stop("input matrix is empty!")
  }
  
  metaData.use = obj@metaData;
  if((x=nrow(metaData.use)) == 0L){
    stop("metaData is empty!")   
  }
  
  ncell = length(obj@barcode)
  nvar = dim(obj@reductions$SVD@x)[2];
  
  
  
  if(any(eigs.dims > nvar) ){
    stop("'eigs.dims' exceeds PCA dimentions number");
  }  
  
  
  
  pca.use = obj@reductions$SVD@x;
  if((x=nrow(pca.use)) == 0L){
    stop("dimentionality reduction is empty, runLDM first")
  }else{
    pca.use = pca.use[,eigs.dims]
  }
  
  
  colnames(x = gmat.use) = paste0(obj@barcode);
  rownames(x = pca.use)  = paste0(obj@barcode);
  rownames(metaData.use) = paste0(obj@barcode);
  seurat.atac <- CreateSeuratObject(counts = data.use, assay = "ATAC");
  seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = as.data.frame(gmat.use))
  seurat.atac <- AddMetaData(seurat.atac, metadata = metaData.use);
  seurat.atac$tech <- "atac"
  DefaultAssay(seurat.atac) <- "ATAC";
  
  colnames(x = pca.use) <- paste0("DC_", eigs.dims);
  seurat.atac[["SnapATAC"]] <- new(Class = "DimReduc", cell.embeddings =pca.use,
                                   feature.loadings = matrix(0,0,0), feature.loadings.projected = matrix(0,0,0),
                                   assay.used ="ATAC", stdev = rep(1,length(eigs.dims)), 
                                   key ="DC_", jackstraw = new(Class = "JackStrawData"), misc = list()) 
  
  DefaultAssay(seurat.atac) <- "ACTIVITY"
  if(norm){ seurat.atac <- NormalizeData(seurat.atac)  }
  if(scale){seurat.atac <- ScaleData(seurat.atac)}
  return(seurat.atac)
}


#' Convert a snap object to a seurat object
#'
#' @param obj A snap object.
#' @param eigs.dims A vector of the dimensions to use.
#' @param norm A logical variable indicates whether to normalize the cell by gene matrix.
#' @param scale A logical variable indicagtes whether to scale the cell by gene matrix.
#' @return Return a seurat object that contains single cell ATAC-seq data
#'
#' @export

art2snap <- function(
  obj, 
  input.mat ="bmat"){
  cat("Epoch: checking input parameters ... \n", file = stderr())
  if(missing(obj)){
    stop("obj is missing")
  }else{
    if(!is(obj,'scART')){
      stop("obj is not a snap object");
    }  
    if((x=length(obj@barcode))==0L){
      stop("obj@barcode is empty");   
    }
  }
  # check if Seurat is installed
  if (requireNamespace("SnapATAC", quietly = TRUE)) {
    require(SnapATAC)
  } else {
    stop("Please install SnapATAC");
  }
  
  if((x=nrow(obj@gmat)) == 0L){
    stop("gmat in obj is empty!");
  }
  gmat.use = (as(obj@gmat,'dgCMatrix'))
  # check if mat is binary;
  
  
  if(input.mat == "bmat"){
    data.use = obj@bmat$filter
    peak.use = as.data.frame(obj@feature);
  }else{
    data.use = obj@gmat;
    peak.use = as.data.frame(rownames(data.use));
  }
  
  if((x=nrow(data.use)) == 0L){
    stop("input matrix is empty!")
  }
  
  metaData.use = obj@metaData;
  if((x=nrow(metaData.use)) == 0L){
    stop("metaData is empty!")   
  }
  
  ncell = length(obj@barcode)
  nvar = dim(obj@reductions$SVD@x)[2];
  
  
  
  if(any(eigs.dims > nvar) ){
    stop("'eigs.dims' exceeds PCA dimentions number");
  }  
  
  
  
  pca.use = obj@reductions$SVD@x;
  if((x=nrow(pca.use)) == 0L){
    stop("dimentionality reduction is empty, runLDM first")
  }else{
    pca.use = pca.use[,eigs.dims]
  }
  
  
  colnames(x = gmat.use) = paste0(obj@barcode);
  rownames(x = pca.use)  = paste0(obj@barcode);
  rownames(metaData.use) = paste0(obj@barcode);
  
  
  snap<- createSnapFromBmat(t(data.use),barcodes = obj@barcode,bins = obj@feature)
  snap@gmat<-gmat.use
  # if(norm){ pbmc.atac <- NormalizeData(pbmc.atac)  }
  # if(scale){pbmc.atac <- ScaleData(pbmc.atac)} 
  return(snap)
}


