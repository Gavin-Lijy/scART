library(readxl)
setwd('./GRN')
regulon<-read.csv( "adj_add_cor_new.csv")
regulon<-regulon[regulon$regulation==1,]
TF.use=unique(regulon$TF)



regulon<-regulon %>% group_by(TF)%>% top_n(20)

library(MotifDb)
library(universalmotif)
library(TFBSTools)
library(motifmatchr)
library(parallel)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg19)
genome=BSgenome.Hsapiens.UCSC.hg19
module=c()
suppressWarnings(annotation <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86))
seqlevelsStyle(annotation) <- "UCSC"
links<-ComputeDORC(peak.data = t(art@pmat),expression.data=art@gmat,annotation = annotation,genome = genome,genes.use = unique(regulon$target))
module=pblapply(1:length(TF.use), function(i){
TF=TF.use[i]
target<-regulon[regulon$TF==TF,]$target
motifs <- query(MotifDb, andStrings=c(TF, "Hsapien"))
a=convert_motifs(motifs, "TFBSTools-PFMatrix")
m=do.call(PFMatrixList,a)
links_TF=links[links$gene %in% target,]
motif_ix=matchMotifs(m, links_TF, genome = genome)
links_TF=links_TF[ !(rowSums( motif_ix@assays@data$motifMatches)==0),]
module=unfactor(unique(links_TF$gene))
return(module)
})
module=lapply(1:length(module), function(i){unfactor(unlist(module[i]))})
exprMat=art@gmat
nCores=1
runSCENIC_3_scoreCells <- function(scenicOptions, exprMat, 
                                   skipBinaryThresholds=FALSE, skipHeatmap=FALSE, skipTsne=FALSE)

  
names(module)=TF.use
regulons=module
  ################################################################
  ## Prepare regulons
 
  regulons <- regulons[order(lengths(regulons), decreasing=TRUE)]
  
  if(length(regulons) <2)  stop("Not enough regulons with at least 10 genes.")
  
  # Add the TF to the regulon (keeping it only once) & rename regulon
  regulons <- setNames(lapply(names(regulons), function(tf) sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))), names(regulons))
  # names(regulons) <- paste(names(regulons), " (",lengths(regulons), "g)", sep="")
  # saveRDS(regulons, file=getIntName(scenicOptions, "aucell_regulons"))
  
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 3. Analyzing the network activity in each individual cell")
  message(msg)
  
  biggestRegulons <- grep("_extended",names(regulons),invert = T, value = T)
  biggestRegulons <- biggestRegulons[1:min(length(biggestRegulons),10)]
  msg <- paste0("\tNumber of regulons to evaluate on cells: ", length(regulons),
                "\nBiggest (non-extended) regulons: \n",
                paste("\t", biggestRegulons, collapse="\n")) # TODO maxlen?
  message(msg)
  
  ################################################################
  # AUCell
  library(AUCell)
  # 1. Create rankings
 
 
  aucellRankings <- AUCell_buildRankings(exprMat, nCores=nCores, 
                                           plotStats=TRUE)
    

  # saveRDS(aucellRankings, file=getIntName(scenicOptions, "aucell_rankings"))
  
  # 2. Calculate AUC
  regulonAUC <- AUCell_calcAUC(regulons, aucellRankings, 
                              nCores=nCores)
  
  # # Order the modules by similarity, for easier exploration in the upcoming steps & save
  # variableRegulons <- names(which(apply(getAUC(regulonAUC), 1, sd) > 0))
  # reguDist <-as.dist(1-cor(t(getAUC(regulonAUC)[variableRegulons,]), method="spear"))
  # reguClust <- hclust(reguDist, method="ward.D2")
  # regulonClusters <- setNames(dynamicTreeCut::cutreeDynamic(reguClust, distM=as.matrix(reguDist), verbose = FALSE), reguClust$labels)
  # regulonOrder <- reguClust$labels[reguClust$order]
  # regulonOrder <- regulonOrder[order(regulonClusters[regulonOrder], decreasing = TRUE)]
  # 
  # regulonAUC <- regulonAUC[regulonOrder,]
  auc= as(regulonAUC@assays@data@listData[["AUC"]],'Matrix')
  art@rmat=auc

  