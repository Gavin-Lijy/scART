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
#' @name scART-class
#' @rdname scART-class
#' @exportClass scART
#' @importFrom methods setClassUnion GenomicRanges Matrix
#' @import   GenomicRanges Matrix
library(methods)
library(GenomicRanges)
setClass('scART',slots=list(barcode="character",feature='GRanges',metaData="data.frame",bmat = "list",smat='Matrix',gmat = "Matrix",
                            mmat = "Matrix",reductions = "list",trajectory='Matrix' ))

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
setMethod("show", signature = "scART",
          definition = function(object) {
            cat("number of barcodes: ", ifelse(is.null(length(object@barcode)), 0, length(object@barcode)), "\n", sep="");
            cat("number of bins: ", ncol(object@bmat), "\n", sep="");
            cat("number of genes: ", ncol(object@gmat), "\n", sep="");
            cat("number of motifs: ", ncol(object@mmat), "\n", sep="");
          }
)
setClass( 'SVD',slots = c(x='matrix',sdev='numeric'))
setClass('TSNE',slots = c(matrix='matrix',nSV='numeric'))
setClass('TSNE_3D',slots = c(matrix='matrix',nSV='numeric'))
setMethod("show", signature = "scART",
          definition = function(object) {
            cat("number of barcodes: ", ifelse(is.null(length(object@barcode)), 0, length(object@barcode)), "\n", sep="");
            cat("number of bins: ", ncol(object@bmat), "\n", sep="");
            cat("number of genes: ", ncol(object@gmat), "\n", sep="");
            cat("number of motifs: ", ncol(object@mmat), "\n", sep="");
          }
)
