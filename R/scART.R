#' An S4 class scART to represent single-nucleus accessibility object.
#'
#' Class art defines a scART object.
#'
#' @slot barcode A character vector contains cell barcodes as rows.
#' @slot feature A matrix object contains the bed file.
#' @slot metaData A data.frame object contains meta data for barcodes.
#' @slot bmat A Matrix object contains all cellxbin count matrixs.
#' @slot smat A matrix object contains the similarity matrix
#' @slot reductions A list object contains the all the reduction results.
#' @slot mmat A matrix object contains cellxmotif variability matrix.
#' @slot trajectory A matrix object contains DDRTREE for trajectory.
#' @slot gmat A matrix object contains cellxgene count matrix.
#' @slot mmat A matrix object contains cellxmotif variability matrix.
#' @slot pmat A Matrix object contains cellxpeak count matrix.
#' @name scART-class
#' @rdname scART-class
#' @exportClass scART
#' @importFrom methods setClassUnion 
#' @import  Matrix GenomicRanges
library(methods)
library(Matrix)
library(GenomicRanges)
setClass('scART',slots=list(barcode="character",feature='GRanges',metaData="data.frame",
                            bmat = "list",smat='Matrix',gmat = "Matrix",
                            mmat = "Matrix", pmat="Matrix",reductions = "list",trajectory='Matrix' ))#,peak='GRanges'

.valid.scART.barcode <- function(object)
{
  if(length(object@barcode) != nrow(object@metaData)){
    return("slot 'barcode' have different length from 'metaData'")
  }
  NULL
}

.valid.scART <- function(object)
{
  
  c(.valid.scART.barcode(object))
}
# methods::setValidity("scART", .valid.scART)

setClass( 'SVD',slots = c(x='matrix',sdev='numeric'))
setClass('TSNE',slots = c(matrix='matrix',nSV='numeric'))
setClass('TSNE_3D',slots = c(matrix='matrix',nSV='numeric'))
setMethod("show", signature = "scART",
          definition = function(object) {
            cat("number of barcodes: ", ifelse(is.null(length(object@barcode)), 0, length(object@barcode)), "\n", sep="");
            cat("number of original bins: ", nrow(object@bmat$raw), "\n", sep="");
            if(!(is.null(object@bmat$filter))){cat("number of filtered bins: ", nrow(object@bmat$filter), "\n", sep="");}
            cat("number of genes: ", nrow(object@gmat), "\n", sep="");
            cat("number of motifs: ", nrow(object@mmat), "\n", sep="");
          }
)
setMethod("[", "scART",
          function(x,i,j,mat=c("bmat", "pmat", "gmat"), drop="missing"){
            .barcode = x@barcode;
            .feature = x@feature;
            .bmat = x@bmat;   
            .gmat = x@gmat;
            .mmat = x@mmat;
            .smat = x@smat;
            .pmat = x@pmat;
            .reductions = x@reductions;
            .trajectory = x@trajectory
            .metaData = x@metaData;
            # a single row or column
            if(!missing(i)){
              if(max(i) > nrow(x)){
                stop("idx exceeds number of cells");
              }
              if(nrow(.bmat) > 0){.bmat <- .bmat[i,,drop=FALSE]}
              if(nrow(.pmat) > 0){.pmat <- .pmat[i,,drop=FALSE]}
              if(nrow(.gmat) > 0){.gmat <- .gmat[i,,drop=FALSE]}	   
              if(nrow(.mmat) > 0){.mmat <- .mmat[i,,drop=FALSE]}	   
              if(nrow(.jmat@jmat) > 0){.jmat <- .jmat[i,,drop=FALSE]}
              if(nrow(.smat@dmat) > 0){.smat <- .smat[i,,drop=FALSE]}
              if(nrow(.tsne) > 0){.tsne <- .tsne[i,,drop=FALSE]}
              if(nrow(.umap) > 0){.umap <- .umap[i,,drop=FALSE]}
              if(nrow(.graph@mat) > 0){.graph <- .graph[i,,drop=FALSE]}
              if(nrow(.metaData) > 0){.metaData <- .metaData[i,,drop=FALSE]}
              if(length(.cluster) > 0){.cluster <- .cluster[i,drop=FALSE]}
              if(length(.barcode) > 0){.barcode <- .barcode[i,drop=FALSE]}
              if(length(.file) > 0){.file <- .file[i,drop=FALSE]}
              if(length(.sample) > 0){.sample <- .sample[i,drop=FALSE]}
            }
            if(!missing(j)){
              mat = match.arg(mat);
              if(mat == "bmat"){
                if(ncol(.bmat) > 0){.bmat <- .bmat[,j,drop=FALSE]}
                if(length(.feature) > 0){.feature <- .feature[j];}	   
              }else if(mat == "pmat"){
                if(ncol(.pmat) > 0){.pmat <- .pmat[,j,drop=FALSE]}
                if(length(.peak) > 0){.peak <- .peak[j];}	   
              }else if(mat == "gmat"){
                if(ncol(.gmat) > 0){.gmat <- .gmat[,j,drop=FALSE]}
              }
            }
            x@bmat = .bmat;
            x@pmat = .pmat;
            x@gmat = .gmat;
            x@mmat = .mmat;
            x@barcode = .barcode;
            x@file = .file;
            x@sample = .sample;
            x@peak = .peak;
            x@feature = .feature;
            x@metaData = .metaData;
            x@umap = .umap;
            x@feature = .feature;
            x@jmat = .jmat;
            x@smat = .smat;
            x@graph = .graph;
            x@cluster = .cluster;
            x@tsne = .tsne;
            return(x);
          })
