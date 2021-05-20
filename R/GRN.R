
gmat=FilterGmat<-function(obj,hvg=2000){
  gmat=obj@gmat[!is.na(rownames(gmat)),]
  library(Seurat)
  seu=CreateSeuratObject(counts = gmat)
  seu=NormalizeData(seu)
  seu=FindVariableFeatures(seu,nfeatures = hvg)
  VariableFeaturePlot(seu)
  gmat=seu@assays$RNA@data
  return(gmat)
}

runGRNboost=function(gmat,path.to.pyscenic,gmat_file='gmat.csv',tf_file,adj_file='adj.csv',out_file='adj_add_cor.csv',min_importance=2,only_pos=TRUE){
  write.csv(gmat,file = gmat_file)
  system(paste(path.to.pyscenic,'grn',gmat_file,tf_file,'-o',adj_file))
  system(paste(path.to.pyscenic,'add_cor',adj_file,gmat_file,'-o',out_file))
  regulon<-read.csv(out_file)
  if(only_pos==TRUE){
    regulon<-regulon[regulon$regulation==1,]
  }
  regulon<-regulon[regulon$importance>min_importance,]
  return(regulon)
}

