
Read_10X <- function(
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
    peaks.loc <- file.path(run, 'peaks.bed')
    matrix.loc <- file.path(run, 'matrix.mtx')
    # Flag to indicate if this data is from CellRanger >= 3.0
     
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!file.exists(peaks.loc) ) {
      stop("Peaks name or features file missing. Expecting ", basename(path = features.loc))
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

  
Read_10X_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- names(x = infile)
  output <- list()
  if (hdf5r::existsGroup(infile, 'matrix')) {
    # cellranger version 3
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      giveCsparse = FALSE
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  output$matrix@Dimnames[[1]]=gsub(':','-', output$matrix@Dimnames[[1]])
  obj<-CreatescART(output$matrix)
    return(obj)
  }

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

Read_snap<-function(file,sample="atac"){
  library(SnapATAC)
  obj<-createSnap(file = file,sample = sample)
  data<-demo.sp@bmat
  data@Dimnames[[1]] <- gsub(':','-',obj@feature$name) 
  data@Dimnames[[2]]<- obj@barcode
  art<-CreatescART(data)
  return(art)
}

