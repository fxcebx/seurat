#' @include generics.R
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix colSums rowSums colMeans rowMeans
#' @importFrom methods setClass setOldClass setClassUnion slot
#' slot<- setMethod new signature slotNames is
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib Seurat
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setOldClass(Classes = 'package_version')
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'OptionalCharacter', members = c('NULL', 'character'))

#' The AnchorSet Class
#'
#' The AnchorSet class is an intermediate data storage class that stores the anchors and other
#' related information needed for performing downstream analyses - namely data integration
#' (\code{\link{IntegrateData}}) and data transfer (\code{\link{TransferData}}).
#'
#' @slot object.list List of objects used to create anchors
#' @slot reference.cells List of cell names in the reference dataset - needed when performing data
#' transfer.
#' @slot reference.objects Position of reference object/s in object.list
#' @slot query.cells List of cell names in the query dataset - needed when performing data transfer
#' @slot anchors The anchor matrix. This contains the cell indices of both anchor pair cells, the
#' anchor score, and the index of the original dataset in the object.list for cell1 and cell2 of
#' the anchor.
#' @slot offsets The offsets used to enable cell look up in downstream functions
#' @slot anchor.features The features used when performing anchor finding.
#' @slot command Store log of parameters that were used
#'
#' @name AnchorSet-class
#' @rdname AnchorSet-class
#' @exportClass AnchorSet
#'
AnchorSet <- setClass(
  Class = "AnchorSet",
  slots = list(
    object.list = "list",
    reference.cells = "vector",
    reference.objects = "vector",
    query.cells = "vector",
    anchors = "ANY",
    offsets = "ANY",
    anchor.features = "ANY",
    command = "ANY"
  )
)

#' The Assay Class
#'
#' The Assay object is the basic unit of Seurat; each Assay stores raw, normalized, and scaled data
#' as well as cluster information, variable features, and any other assay-specific metadata.
#' Assays should contain single cell expression data such as RNA-seq, protein, or imputed expression
#' data.
#'
#' @slot counts Unnormalized data such as raw counts or TPMs
#' @slot data Normalized expression data
#' @slot scale.data Scaled expression data
#' @slot key Key for the Assay
#' @slot assay.orig Original assay that this assay is based off of. Used to track
#' assay provenence
#' @slot var.features Vector of features exhibiting high variance across single cells
#' @slot meta.features Feature-level metadata
#' @slot misc Utility slot for storing additional data associated with the assay
#'
#' @name Assay-class
#' @rdname Assay-class
#' @exportClass Assay
#'
Assay <- setClass(
  Class = 'Assay',
  slots = c(
    counts = 'AnyMatrix',
    data = 'AnyMatrix',
    scale.data = 'matrix',
    key = 'character',
    assay.orig = 'OptionalCharacter',
    var.features = 'vector',
    meta.features = 'data.frame',
    misc = 'ANY'
  )
)

#' The JackStrawData Class
#'
#' The JackStrawData is used to store the results of a JackStraw computation.
#'
#' @slot empirical.p.values Empirical p-values
#' @slot fake.reduction.scores Fake reduction scores
#' @slot empirical.p.values.full Empirical p-values on full
#' @slot overall.p.values Overall p-values from ScoreJackStraw
#'
#' @name JackStrawData-class
#' @rdname JackStrawData-class
#' @exportClass JackStrawData
#'
JackStrawData <- setClass(
  Class = "JackStrawData",
  slots = list(
    empirical.p.values = "matrix",
    fake.reduction.scores = "matrix",
    empirical.p.values.full = "matrix",
    overall.p.values = "matrix"
  )
)

#' The Dimmensional Reduction Class
#'
#' The DimReduc object stores a dimensionality reduction taken out in Seurat; each DimReduc
#' consists of a cell embeddings matrix, a feature loadings matrix, and a projected feature
#' loadings matrix.
#'
#' @slot cell.embeddings Cell embeddings matrix (required)
#' @slot feature.loadings Feature loadings matrix (optional)
#' @slot feature.loadings.projected Projected feature loadings matrix (optional)
#' @slot assay.used Name of assay used to generate \code{DimReduc} object
#' @slot global Is this \code{DimReduc} global/persistent? If so, it will not be
#' removed when removing its associated assay
#' @slot stdev A vector of standard deviations
#' @slot key Key for the \code{DimReduc}, must be alphanumerics followed by an underscore
#' @slot jackstraw A \code{\link{JackStrawData-class}} object associated with
#' this \code{DimReduc}
#' @slot misc Utility slot for storing additional data associated with the
#' \code{DimReduc} (e.g. the total variance of the PCA)
#'
#' @name DimReduc-class
#' @rdname DimReduc-class
#' @exportClass DimReduc
#'
DimReduc <- setClass(
  Class = 'DimReduc',
  slots = c(
    cell.embeddings = 'matrix',
    feature.loadings = 'matrix',
    feature.loadings.projected = 'matrix',
    assay.used = 'character',
    global = 'logical',
    stdev = 'numeric',
    key = 'character',
    jackstraw = 'JackStrawData',
    misc = 'list'
  )
)

#' The Graph Class
#'
#' The Graph class inherits from dgCMatrix. We do this to enable future expandability of graphs.
#'
#' @slot assay.used Optional name of assay used to generate \code{Graph} object
#'
#' @name Graph-class
#' @rdname Graph-class
#' @exportClass Graph
#'
#' @seealso \code{\link[Matrix]{dgCMatrix-class}}
#'
Graph <- setClass(
  Class = 'Graph',
  contains = "dgCMatrix",
  slots = list(
    assay.used = 'OptionalCharacter'
  )
)
#' The SeuratCommand Class
#'
#' The SeuratCommand is used for logging commands that are run on a SeuratObject. It stores parameters and timestamps
#'
#' @slot name Command name
#' @slot time.stamp Timestamp of when command was tun
#' @slot assay.used Optional name of assay used to generate \code{SeuratCommand} object
#' @slot call.string String of the command call
#' @slot params List of parameters used in the command call
#'
#' @name SeuratCommand-class
#' @rdname SeuratCommand-class
#' @exportClass SeuratCommand
#'
SeuratCommand <- setClass(
  Class = 'SeuratCommand',
  slots = c(
    name = 'character',
    time.stamp = 'POSIXct',
    assay.used = 'OptionalCharacter',
    call.string = 'character',
    params = 'ANY'
  )
)

#' The Seurat Class
#'
#' The Seurat object is a representation of single-cell expression data for R; each Seurat
#' object revolves around a set of cells and consists of one or more \code{\link{Assay-class}}
#' objects, or individual representations of expression data (eg. RNA-seq, ATAC-seq, etc).
#' These assays can be reduced from their high-dimensional state to a lower-dimension state
#' and stored as \code{\link{DimReduc-class}} objects. Seurat objects also store additional
#' meta data, both at the cell and feature level (contained within individual assays). The
#' object was designed to be as self-contained as possible, and easily extendible to new methods.
#'
#' @slot assays A list of assays for this project
#' @slot meta.data Contains meta-information about each cell, starting with number of genes detected (nGene)
#' and the original identity class (orig.ident); more information is added using \code{AddMetaData}
#' @slot active.assay Name of the active, or default, assay; settable using \code{\link{DefaultAssay}}
#' @slot active.ident The active cluster identity for this Seurat object; settable using \code{\link{Idents}}
#' @slot graphs A list of \code{\link{Graph-class}} objects
#' @slot neighbors ...
#' @slot reductions A list of dimmensional reduction objects for this object
#' @slot project.name Name of the project
#' @slot misc A list of miscellaneous information
#' @slot version Version of Seurat this object was built under
#' @slot commands A list of logged commands run on this \code{Seurat} object
#' @slot tools A list of miscellaneous data generated by other tools, should be filled by developers only using \code{\link{Tool}<-}
#'
#' @name Seurat-class
#' @rdname Seurat-class
#' @exportClass Seurat
#'
Seurat <- setClass(
  Class = 'Seurat',
  slots = c(
    assays = 'list',
    meta.data = 'data.frame',
    active.assay = 'character',
    active.ident = 'factor',
    graphs = 'list',
    neighbors = 'list',
    reductions = 'list',
    project.name = 'character',
    misc = 'list',
    version = 'package_version',
    commands = 'list',
    tools = 'list'
  )
)

#' The Seurat Class
#'
#' The Seurat object is the center of each single cell analysis. It stores all information
#' associated with the dataset, including data, annotations, analyes, etc. All that is needed
#' to construct a Seurat object is an expression matrix (rows are genes, columns are cells), which
#' should be log-scale
#'
#' Each Seurat object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot raw.data The raw project data
#' @slot data The normalized expression matrix (log-scale)
#' @slot scale.data scaled (default is z-scoring each gene) expression matrix; used for dimmensional reduction and heatmap visualization
#' @slot var.genes Vector of genes exhibiting high variance across single cells
#' @slot is.expr Expression threshold to determine if a gene is expressed (0 by default)
#' @slot ident THe 'identity class' for each cell
#' @slot meta.data Contains meta-information about each cell, starting with number of genes detected (nGene)
#' and the original identity class (orig.ident); more information is added using \code{AddMetaData}
#' @slot project.name Name of the project (for record keeping)
#' @slot dr List of stored dimmensional reductions; named by technique
#' @slot assay List of additional assays for multimodal analysis; named by technique
#' @slot hvg.info The output of the mean/variability analysis for all genes
#' @slot imputed Matrix of imputed gene scores
#' @slot cell.names Names of all single cells (column names of the expression matrix)
#' @slot cluster.tree List where the first element is a phylo object containing the phylogenetic tree relating different identity classes
#' @slot snn Spare matrix object representation of the SNN graph
#' @slot calc.params Named list to store all calculation-related parameter choices
#' @slot kmeans Stores output of gene-based clustering from \code{DoKMeans}
#' @slot spatial Stores internal data and calculations for spatial mapping of single cells
#' @slot misc Miscellaneous spot to store any data alongisde the object (for example, gene lists)
#' @slot version Version of package used in object creation
#'
#' @name seurat-class
#' @rdname oldseurat-class
#' @aliases seurat-class
#'
seurat <- setClass(
  Class = "seurat",
  slots = c(
    raw.data = "ANY",
    data = "ANY",
    scale.data = "ANY",
    var.genes = "vector",
    is.expr = "numeric",
    ident = "factor",
    meta.data = "data.frame",
    project.name = "character",
    dr = "list",
    assay = "list",
    hvg.info = "data.frame",
    imputed = "data.frame",
    cell.names = "vector",
    cluster.tree = "list",
    snn = "dgCMatrix",
    calc.params = "list",
    kmeans = "ANY",
    spatial = "ANY",
    misc = "ANY",
    version = "ANY"
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Pull Assays or assay names
#'
#' Lists the names of \code{\link{Assay}} objects present in
#' a Seurat object. If slot is provided, pulls specified Assay object.
#'
#' @param object A Seurat object
#' @param slot Name of Assay to return
#'
#' @return If \code{slot} is \code{NULL}, the names of all \code{Assay} objects
#' in this Seurat object. Otherwise, the \code{Assay} object specified
#'
#' @export
#'
#' @examples
#' Assays(object = pbmc_small)
#'
Assays <- function(object, slot = NULL) {
  assays <- FilterObjects(object = object, classes.keep = 'Assay')
  if (is.null(x = slot)) {
    return(assays)
  }
  if (!slot %in% assays) {
    warning(
      "Cannot find an assay of name ",
      slot,
      " in this Seurat object",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(slot(object = object, name = 'assays')[[slot]])
}

#' Get cell names grouped by identity class
#'
#' @param object A Seurat object
#' @param idents A vector of identity class levels to limit resulting list to;
#' defaults to all identity class levels
#' @param cells A vector of cells to grouping to
#'
#' @return A named list where names are identity classes and values are vectors
#' of cells beloning to that class
#'
#' @export
#'
#' @examples
#' CellsByIdentities(object = pbmc_small)
#'
CellsByIdentities <- function(object, idents = NULL, cells = NULL) {
  cells <- cells %||% colnames(x = object)
  cells <- intersect(x = cells, y = colnames(x = object))
  if (length(x = cells) == 0) {
    stop("Cannot find cells provided")
  }
  idents <- idents %||% levels(x = object)
  idents <- intersect(x = idents, y = levels(x = object))
  if (length(x = idents) == 0) {
    stop("None of the provided identity class levels were found", call. = FALSE)
  }
  cells.idents <- sapply(
    X = idents,
    FUN = function(i) {
      return(cells[as.vector(x = Idents(object = object)[cells]) == i])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (any(is.na(x = Idents(object = object)[cells]))) {
    cells.idents["NA"] <- names(x = which(x = is.na(x = Idents(object = object)[cells])))
  }
  return(cells.idents)
}

#' Create an Assay object
#'
#' Create an Assay object from a feature (e.g. gene) expression matrix. The
#' expected format of the input matrix is features x cells.
#'
#' Non-unique cell or feature names are not allowed. Please make unique before
#' calling this function.
#'
#' @param counts Unnormalized data such as raw counts or TPMs
#' @param data Prenormalized data; if provided, do not pass \code{counts}
#' @param min.cells Include features detected in at least this many cells. Will
#' subset the counts matrix as well. To reintroduce excluded features, create a
#' new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are
#' detected.
#'
#' @importFrom methods as
#' @importFrom Matrix colSums rowSums
#'
#' @export
#'
#' @examples
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_rna <- CreateAssayObject(counts = pbmc_raw)
#' pbmc_rna
#'
CreateAssayObject <- function(
  counts,
  data,
  min.cells = 0,
  min.features = 0
) {
  if (missing(x = counts) && missing(x = data)) {
    stop("Must provide either 'counts' or 'data'")
  } else if (!missing(x = counts) && !missing(x = data)) {
    stop("Either 'counts' or 'data' must be missing; both cannot be provided")
  } else if (!missing(x = counts)) {
    # check that dimnames of input counts are unique
    if (anyDuplicated(rownames(x = counts))) {
      warning(
        "Non-unique features (rownames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = counts) <- make.unique(names = rownames(x = counts))
    }
    if (anyDuplicated(colnames(x = counts))) {
      warning(
        "Non-unique cell names (colnames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      colnames(x = counts) <- make.unique(names = colnames(x = counts))
    }
    if (is.null(x = colnames(x = counts))) {
      stop("No cell names (colnames) names present in the input matrix")
    }
    if (any(rownames(x = counts) == '')) {
      stop("Feature names of counts matrix cannot be empty", call. = FALSE)
    }
    if (nrow(x = counts) > 0 && is.null(x = rownames(x = counts))) {
      stop("No feature names (rownames) names present in the input matrix")
    }
    if (!inherits(x = counts, what = 'dgCMatrix')) {
      counts <- as(object = as.matrix(x = counts), Class = 'dgCMatrix')
    }
    # Filter based on min.features
    if (min.features > 0) {
      nfeatures <- Matrix::colSums(x = counts > 0)
      counts <- counts[, which(x = nfeatures >= min.features)]
    }
    # filter genes on the number of cells expressing
    if (min.cells > 0) {
      num.cells <- Matrix::rowSums(x = counts > 0)
      counts <- counts[which(x = num.cells >= min.cells), ]
    }
    data <- counts
  } else if (!missing(x = data)) {
    # check that dimnames of input data are unique
    if (anyDuplicated(rownames(x = data))) {
      warning(
        "Non-unique features (rownames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = data) <- make.unique(names = rownames(x = data))
    }
    if (anyDuplicated(colnames(x = data))) {
      warning(
        "Non-unique cell names (colnames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      colnames(x = data) <- make.unique(names = colnames(x = data))
    }
    if (is.null(x = colnames(x = data))) {
      stop("No cell names (colnames) names present in the input matrix")
    }
    if (any(rownames(x = data) == '')) {
      stop("Feature names of data matrix cannot be empty", call. = FALSE)
    }
    if (nrow(x = data) > 0 && is.null(x = rownames(x = data))) {
      stop("No feature names (rownames) names present in the input matrix")
    }
    if (min.cells != 0 | min.features != 0) {
      warning(
        "No filtering performed if passing to data rather than counts",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    counts <- new(Class = 'matrix')
  }
  # Ensure row- and column-names are vectors, not arrays
  if (!is.vector(x = rownames(x = counts))) {
    rownames(x = counts) <- as.vector(x = rownames(x = counts))
  }
  if (!is.vector(x = colnames(x = counts))) {
    colnames(x = counts) <- as.vector(x = colnames(x = counts))
  }
  if (!is.vector(x = rownames(x = data))) {
    rownames(x = data) <- as.vector(x = rownames(x = data))
  }
  if (!is.vector(x = colnames(x = data))) {
    colnames(x = data) <- as.vector(x = colnames(x = data))
  }
  if (any(grepl(pattern = '_', x = rownames(x = counts))) || any(grepl(pattern = '_', x = rownames(x = data)))) {
    warning(
      "Feature names cannot have underscores ('_'), replacing with dashes ('-')",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(x = counts) <- gsub(
      pattern = '_',
      replacement = '-',
      x = rownames(x = counts)
    )
    rownames(x = data) <- gsub(
      pattern = '_',
      replacement = '-',
      x = rownames(x = data)
    )
  }
  if (any(grepl(pattern = '|', x = rownames(x = counts), fixed = TRUE)) || any(grepl(pattern = '|', x = rownames(x = data), fixed = TRUE))) {
    warning(
      "Feature names cannot have pipe characters ('|'), replacing with dashes ('-')",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(x = counts) <- gsub(
      pattern = '|',
      replacement = '-',
      x = rownames(x = counts),
      fixed = TRUE
    )
    rownames(x = data) <- gsub(
      pattern = '|',
      replacement = '-',
      x = rownames(x = data),
      fixed = TRUE
    )
  }
  # Initialize meta.features
  init.meta.features <- data.frame(row.names = rownames(x = data))
  assay <- new(
    Class = 'Assay',
    counts = counts,
    data = data,
    scale.data = new(Class = 'matrix'),
    meta.features = init.meta.features
  )
  return(assay)
}

#' Create a DimReduc object
#'
#' @param embeddings A matrix with the cell embeddings
#' @param loadings A matrix with the feature loadings
#' @param projected A matrix with the projected feature loadings
#' @param assay Assay used to calculate this dimensional reduction
#' @param stdev Standard deviation (if applicable) for the dimensional reduction
#' @param key A character string to facilitate looking up features from a
#' specific DimReduc
#' @param global Specify this as a global reduction (useful for visualizations)
#' @param jackstraw Results from the JackStraw function
#' @param misc list for the user to store any additional information associated
#' with the dimensional reduction
#'
#' @aliases SetDimReduction
#'
#' @export
#'
#' @examples
#' data <- GetAssayData(pbmc_small[["RNA"]], slot = "scale.data")
#' pcs <- prcomp(x = data)
#' pca.dr <- CreateDimReducObject(
#'   embeddings = pcs$rotation,
#'   loadings = pcs$x,
#'   stdev = pcs$sdev,
#'   key = "PC",
#'   assay = "RNA"
#' )
#'
CreateDimReducObject <- function(
  embeddings = new(Class = 'matrix'),
  loadings = new(Class = 'matrix'),
  projected = new(Class = 'matrix'),
  assay = NULL,
  stdev = numeric(),
  key = NULL,
  global = FALSE,
  jackstraw = NULL,
  misc = list()
) {
  if (is.null(x = assay)) {
    warning(
      "No assay specified, setting assay as RNA by default.",
      call. = FALSE,
      immediate. = TRUE
    )
    assay <- "RNA"
  }
  # Try to infer key from column names
  if (is.null(x = key) && is.null(x = colnames(x = embeddings))) {
    stop("Please specify a key for the DimReduc object")
  } else if (is.null(x = key)) {
    key <- regmatches(
      x = colnames(x = embeddings),
      m = regexec(pattern = '^[[:alnum:]]+_', text = colnames(x = embeddings))
    )
    key <- unique(x = unlist(x = key, use.names = FALSE))
  }
  if (length(x = key) != 1) {
    stop("Please specify a key for the DimReduc object")
  } else if (!grepl(pattern = '^[[:alnum:]]+_$', x = key)) {
    # New SetKey function
    key <- regmatches(
      x = key,
      m = regexec(pattern = '[[:alnum:]]+', text = key)
    )
    key <- paste0(paste(key, collapse = ''), '_')
    warning(
      "All keys should be one or more alphanumeric characters followed by an underscore '_', setting key to ",
      key,
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # ensure colnames of the embeddings are the key followed by a numeric
  if (is.null(x = colnames(x = embeddings))) {
    warning(
      "No columnames present in cell embeddings, setting to '",
      key,
      "1:",
      ncol(x = embeddings),
      "'",
      call. = FALSE,
      immediate. = TRUE
    )
    colnames(x = embeddings) <- paste0(key, 1:ncol(x = embeddings))
  } else if (!all(grepl(pattern = paste0('^', key, "[[:digit:]]+$"), x = colnames(x = embeddings)))) {
    digits <- unlist(x = regmatches(
      x = colnames(x = embeddings),
      m = regexec(pattern = '[[:digit:]]+$', text = colnames(x = embeddings))
    ))
    if (length(x = digits) != ncol(x = embeddings)) {
      stop("Please ensure all column names in the embeddings matrix are the key plus a digit representing a dimension number")
    }
    colnames(x = embeddings) <- paste0(key, digits)
  }
  if (!IsMatrixEmpty(x = loadings)) {
    if (any(rownames(x = loadings) == '')) {
      stop("Feature names of loadings matrix cannot be empty", call. = FALSE)
    }
    colnames(x = loadings) <- colnames(x = embeddings)
  }
  if (!IsMatrixEmpty(x = projected)) {
    if (any(rownames(x = loadings) == '')) {
      stop("Feature names of projected loadings matrix cannot be empty", call. = FALSE)
    }
    colnames(x = projected) <- colnames(x = embeddings)
  }
  jackstraw <- jackstraw %||% new(Class = 'JackStrawData')
  dim.reduc <- new(
    Class = 'DimReduc',
    cell.embeddings = embeddings,
    feature.loadings = loadings,
    feature.loadings.projected = projected,
    assay.used = assay,
    global = global,
    stdev = stdev,
    key = key,
    jackstraw = jackstraw,
    misc = misc
  )
  return(dim.reduc)
}

#' Create a Seurat object
#'
#' Create a Seurat object from a feature (e.g. gene) expression matrix. The expected format of the
#' input matrix is features x cells.
#'
#'
#' Note: In previous versions (<3.0), this function also accepted a parameter to set the expression
#' threshold for a 'detected' feature (gene). This functionality has been removed to simplify the
#' initialization process/assumptions. If you would still like to impose this threshold for your
#' particular dataset, simply filter the input expression matrix before calling this function.
#'
#' @inheritParams CreateAssayObject
#' @param project Sets the project name for the Seurat object.
#' @param assay Name of the assay corresponding to the initial input data.
#' @param names.field For the initial identity class for each cell, choose this field from the
#' cell's name. E.g. If your cells are named as BARCODE_CLUSTER_CELLTYPE in the input matrix, set
#' names.field to 3 to set the initial identities to CELLTYPE.
#' @param names.delim For the initial identity class for each cell, choose this delimiter from the
#' cell's column name. E.g. If your cells are named as BARCODE-CLUSTER-CELLTYPE, set this to "-" to
#' separate the cell name into its component parts for picking the relevant field.
#' @param meta.data Additional cell-level metadata to add to the Seurat object. Should be a data
#' frame where the rows are cell names and the columns are additional metadata fields.
#'
#' @importFrom utils packageVersion
#' @importFrom Matrix colSums
#' @export
#'
#' @examples
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_small <- CreateSeuratObject(counts = pbmc_raw)
#' pbmc_small
#'
CreateSeuratObject <- function(
  counts,
  project = 'SeuratProject',
  assay = 'RNA',
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
) {
  if (!is.null(x = meta.data)) {
    if (is.null(x = rownames(x = meta.data))) {
      stop("Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
    }
    if (length(x = setdiff(x = rownames(x = meta.data), y = colnames(x = counts)))) {
      warning("Some cells in meta.data not present in provided counts matrix.")
      meta.data <- meta.data[intersect(x = rownames(x = meta.data), y = colnames(x = counts)), ]
    }
    if (is.data.frame(x = meta.data)) {
      new.meta.data <- data.frame(row.names = colnames(x = counts))
      for (ii in 1:ncol(x = meta.data)) {
        new.meta.data[rownames(x = meta.data), colnames(x = meta.data)[ii]] <- meta.data[, ii, drop = FALSE]
      }
      meta.data <- new.meta.data
    }
  }
  assay.data <- CreateAssayObject(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features
  )
  Key(object = assay.data) <- paste0(tolower(x = assay), '_')
  assay.list <- list(assay.data)
  names(x = assay.list) <- assay
  init.meta.data <- data.frame(row.names = colnames(x = assay.list[[assay]]))
  # Set idents
  idents <- factor(x = unlist(x = lapply(
    X = colnames(x = assay.data),
    FUN = ExtractField,
    field = names.field,
    delim = names.delim
  )))
  if (any(is.na(x = idents))) {
    warning("Input parameters result in NA values for initial cell identities. Setting all initial idents to the project name")
  }
  # if there are more than 100 idents, set all idents to ... name
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == length(x = idents)) {
    idents <- rep.int(x = factor(x = project), times = ncol(x = assay.data))
  }
  names(x = idents) <- colnames(x = assay.data)
  object <- new(
    Class = 'Seurat',
    assays = assay.list,
    meta.data = init.meta.data,
    active.assay = assay,
    active.ident = idents,
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  object[['orig.ident']] <- idents
  # Calculate nCount and nFeature
  n.calc <- CalcN(object = assay.data)
  if (!is.null(x = n.calc)) {
    names(x = n.calc) <- paste(names(x = n.calc), assay, sep = '_')
    object[[names(x = n.calc)]] <- n.calc
  }
  if (!is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  return(object)
}

#' Access cellular data
#'
#' Retreives data (feature expression, PCA scores, metrics, etc.) for a set
#' of cells in a Seurat object
#'
#' @param object Seurat object
#' @param vars List of all variables to fetch, use keyword 'ident' to pull identity classes
#' @param cells Cells to collect data for (default is all cells)
#' @param slot Slot to pull feature data for
#'
#' @return A data frame with cells as rows and cellular data as columns
#'
#' @export
#'
#' @examples
#' pc1 <- FetchData(object = pbmc_small, vars = 'PC_1')
#' head(x = pc1)
#' head(x = FetchData(object = pbmc_small, vars = c('groups', 'ident')))
#'
FetchData <- function(object, vars, cells = NULL, slot = 'data') {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  # Get a list of all objects to search through and their keys
  objects.use <- FilterObjects(object = object)
  object.keys <- sapply(X = objects.use, FUN = function(i) {return(Key(object[[i]]))})
  # Find all vars that are keyed
  keyed.vars <- lapply(
    X = object.keys,
    FUN = function(key) {
      if (length(x = key) == 0) {
        return(integer(length = 0L))
      }
      return(grep(pattern = paste0('^', key), x = vars))
    }
  )
  keyed.vars <- Filter(f = length, x = keyed.vars)
  data.fetched <- lapply(
    X = names(x = keyed.vars),
    FUN = function(x) {
      vars.use <- vars[keyed.vars[[x]]]
      key.use <- object.keys[x]
      data.return <- if (inherits(x = object[[x]], what = 'DimReduc')) {
        vars.use <- grep(
          pattern = paste0('^', key.use, '[[:digit:]]+$'),
          x = vars.use,
          value = TRUE
        )
        if (length(x = vars.use) > 0) {
          tryCatch(
            expr = object[[x]][[cells, vars.use, drop = FALSE]],
            error = function(...) {
              return(NULL)
            }
          )
        } else {
          NULL
        }
      } else if (inherits(x = object[[x]], what = 'Assay')) {
        vars.use <- gsub(pattern = paste0('^', key.use), replacement = '', x = vars.use)
        data.assay <- GetAssayData(
          object = object,
          slot = slot,
          assay = x
        )
        vars.use <- vars.use[vars.use %in% rownames(x = data.assay)]
        data.vars <- t(x = as.matrix(data.assay[vars.use, cells, drop = FALSE]))
        if (ncol(data.vars) > 0) {
          colnames(x = data.vars) <- paste0(key.use, vars.use)
        }
        data.vars
      }
      data.return <- as.list(x = as.data.frame(x = data.return))
      return(data.return)
    }
  )
  data.fetched <- unlist(x = data.fetched, recursive = FALSE)
  # Pull vars from object metadata
  meta.vars <- vars[vars %in% colnames(x = object[[]]) & ! vars %in% names(x = data.fetched)]
  data.fetched <- c(data.fetched, object[[meta.vars]][cells, , drop = FALSE])
  # Pull vars from the default assay
  default.vars <- vars[vars %in% rownames(x = GetAssayData(object = object, slot = slot)) & ! vars %in% names(x = data.fetched)]
  data.fetched <- c(
    data.fetched,
    tryCatch(
      expr = as.data.frame(x = t(x = as.matrix(x = GetAssayData(
        object = object,
        slot = slot
      )[default.vars, cells, drop = FALSE]))),
      error = function(...) {
        return(NULL)
      }
    )
  )
  # Pull identities
  if ('ident' %in% vars && !'ident' %in% colnames(x = object[[]])) {
    data.fetched[['ident']] <- Idents(object = object)[cells]
  }
  # Try to find ambiguous vars
  fetched <- names(x = data.fetched)
  vars.missing <- setdiff(x = vars, y = fetched)
  if (length(x = vars.missing) > 0) {
    # Search for vars in alternative assays
    vars.alt <- vector(mode = 'list', length = length(x = vars.missing))
    names(x = vars.alt) <- vars.missing
    for (assay in FilterObjects(object = object, classes.keep = 'Assay')) {
      vars.assay <- Filter(
        f = function(x) {
          features.assay <- rownames(x = GetAssayData(
            object = object,
            assay = assay,
            slot = slot
          ))
          return(x %in% features.assay)
        },
        x = vars.missing
      )
      for (var in vars.assay) {
        vars.alt[[var]] <- append(x = vars.alt[[var]], values = assay)
      }
    }
    # Vars found in multiple alternative assays are truly ambiguous, will not pull
    vars.many <- names(x = Filter(
      f = function(x) {
        return(length(x = x) > 1)
      },
      x = vars.alt
    ))
    if (length(x = vars.many) > 0) {
      warning(
        "Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe: ",
        paste(vars.many, collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
    }
    vars.missing <- names(x = Filter(
      f = function(x) {
        return(length(x = x) != 1)
      },
      x = vars.alt
    ))
    # Pull vars found in only one alternative assay
    # Key this var to highlight that it was found in an alternate assay
    vars.alt <- Filter(
      f = function(x) {
        return(length(x = x) == 1)
      },
      x = vars.alt
    )
    for (var in names(x = vars.alt)) {
      assay <- vars.alt[[var]]
      warning(
        'Could not find ',
        var,
        ' in the default search locations, found in ',
        assay,
        ' assay instead',
        immediate. = TRUE,
        call. = FALSE
      )
      keyed.var <- paste0(Key(object = object[[assay]]), var)
      data.fetched[[keyed.var]] <- as.vector(
        x = GetAssayData(object = object, assay = assay, slot = slot)[var, cells]
      )
      vars <- sub(
        pattern = paste0('^', var, '$'),
        replacement = keyed.var,
        x = vars
      )
    }
    fetched <- names(x = data.fetched)
  }
  # Name the vars not found in a warning (or error if no vars found)
  m2 <- if (length(x = vars.missing) > 10) {
    paste0(' (10 out of ', length(x = vars.missing), ' shown)')
  } else {
    ''
  }
  if (length(x = vars.missing) == length(x = vars)) {
    stop(
      "None of the requested variables were found",
      m2,
      ': ',
      paste(head(x = vars.missing, n = 10L), collapse = ', ')
    )
  } else if (length(x = vars.missing) > 0) {
    warning(
      "The following requested variables were not found",
      m2,
      ': ',
      paste(head(x = vars.missing, n = 10L), collapse = ', ')
    )
  }
  # Assembled fetched vars in a dataframe
  data.fetched <- as.data.frame(
    x = data.fetched,
    row.names = cells,
    stringsAsFactors = FALSE
  )
  data.order <- na.omit(object = pmatch(
    x = vars,
    table = fetched
  ))
  if (length(x = data.order) > 1) {
    data.fetched <- data.fetched[, data.order]
  }
  colnames(x = data.fetched) <- vars[vars %in% fetched]
  return(data.fetched)
}

#' Log a command
#'
#' Logs command run, storing the name, timestamp, and argument list. Stores in
#' the Seurat object
#'
#' @param object Name of Seurat object
#' @param return.command Return a \link{SeuratCommand} object instead
#'
#' @return If \code{return.command}, returns a SeuratCommand object. Otherwise,
#' returns the Seurat object with command stored
#'
#' @export
#'
#' @seealso \code{\link{Command}}
#'
LogSeuratCommand <- function(object, return.command = FALSE) {
  time.stamp <- Sys.time()
  #capture function name
  which.frame <- sys.nframe() - 1
  if (which.frame < 1) {
    stop("'LogSeuratCommand' cannot be called at the top level", call. = FALSE)
  }
  command.name <- as.character(x = deparse(expr = sys.calls()[[which.frame]]))
  command.name <- gsub(pattern = "\\.Seurat", replacement = "", x = command.name)
  call.string <- command.name
  command.name <- ExtractField(string = command.name, field = 1, delim = "\\(")
  #capture function arguments
  argnames <- names(x = formals(fun = sys.function(which = sys.parent(n = 1))))
  argnames <- grep(pattern = "object", x = argnames, invert = TRUE, value = TRUE)
  argnames <- grep(pattern = "anchorset", x = argnames, invert = TRUE, value = TRUE)
  argnames <- grep(pattern = "\\.\\.\\.", x = argnames, invert = TRUE, value = TRUE)
  params <- list()
  p.env <- parent.frame(n = 1)
  argnames <- intersect(x = argnames, y = ls(name = p.env))
  # fill in params list
  for (arg in argnames) {
    param_value <- get(x = arg, envir = p.env)
    if (inherits(x = param_value, what = 'Seurat')) {
      next
    }
    #TODO Institute some check of object size?
    params[[arg]] <- param_value
  }
  # check if function works on the Assay and/or the DimReduc Level
  assay <- params[["assay"]]
  reduction <- params[["reduction"]]
  # Get assay used for command
  cmd.assay <- assay %||% (reduction %iff% if (inherits(x = reduction, what = 'DimReduc')) {
    DefaultAssay(object = reduction)
  } else if (reduction %in% Reductions(object = object)) {
    DefaultAssay(object = object[[reduction]])
  })
  if (inherits(x = reduction, what = 'DimReduc')) {
    reduction <- 'DimReduc'
  }
  # rename function name to include Assay/DimReduc info
  if (length(x = assay) == 1) {
    command.name <- paste(command.name, assay, reduction, sep = '.')
  }
  command.name <- sub(pattern = "[\\.]+$", replacement = "", x = command.name, perl = TRUE)
  command.name <- sub(pattern = "\\.\\.", replacement = "\\.", x = command.name, perl = TRUE)
  # store results
  seurat.command <- new(
    Class = 'SeuratCommand',
    name = command.name,
    params = params,
    time.stamp = time.stamp,
    call.string = call.string,
    assay.used = cmd.assay
  )
  if (return.command) {
    return(seurat.command)
  }
  object[[command.name]] <- seurat.command
  return(object)
}

#' Pull DimReducs or DimReduc names
#'
#' Lists the names of \code{\link{DimReduc}} objects present in
#' a Seurat object. If slot is provided, pulls specified DimReduc object.
#'
#' @param object A Seurat object
#' @param slot Name of DimReduc
#'
#' @return If \code{slot} is \code{NULL}, the names of all \code{DimReduc} objects
#' in this Seurat object. Otherwise, the \code{DimReduc} object requested
#'
#' @export
#'
#' @examples
#' Reductions(object = pbmc_small)
#'
Reductions <- function(object, slot = NULL) {
  reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
  if (is.null(x = slot)) {
    return(reductions)
  }
  if (!slot %in% reductions) {
    warning(
      "Cannot find a DimReduc of name ",
      slot,
      " in this Seurat object",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(slot(object = object, name = 'reductions')[[slot]])
}

#' Find features with highest scores for a given dimensional reduction technique
#'
#' Return a list of features with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim Dimension to use
#' @param nfeatures Number of features to return
#' @param projected Use the projected feature loadings
#' @param balanced Return an equal number of features with both + and - scores.
#' @param ... Extra parameters passed to \code{\link{Loadings}}
#'
#' @return Returns a vector of features
#'
#' @export
#'
#' @examples
#' pbmc_small
#' TopFeatures(object = pbmc_small[["pca"]], dim = 1)
#' # After projection:
#' TopFeatures(object = pbmc_small[["pca"]], dim = 1,  projected = TRUE)
#'
TopFeatures <- function(
  object,
  dim = 1,
  nfeatures = 20,
  projected = FALSE,
  balanced = FALSE,
  ...
) {
  loadings <- Loadings(object = object, projected = projected, ...)[, dim, drop = FALSE]
  return(Top(
    data = loadings,
    num = nfeatures,
    balanced = balanced
  ))
}

#' Find cells with highest scores for a given dimensional reduction technique
#'
#' Return a list of genes with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim Dimension to use
#' @param ncells Number of cells to return
#' @param balanced Return an equal number of cells with both + and - scores.
#' @param ... Extra parameters passed to \code{\link{Embeddings}}
#'
#' @return Returns a vector of cells
#'
#' @export
#'
#' @examples
#' pbmc_small
#' head(TopCells(object = pbmc_small[["pca"]]))
#' # Can specify which dimension and how many cells to return
#' TopCells(object = pbmc_small[["pca"]], dim = 2, ncells = 5)
#'
TopCells <- function(object, dim = 1, ncells = 20, balanced = FALSE, ...) {
  embeddings <- Embeddings(object = object, ...)[, dim, drop = FALSE]
  return(Top(
    data = embeddings,
    num = ncells,
    balanced = balanced
  ))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname AddMetaData
#' @export
#' @method AddMetaData Seurat
#'
AddMetaData.Seurat <- function(object, metadata, col.name = NULL) {
  return(.AddMetaData(object = object, metadata = metadata, col.name = col.name))
}

#' @param assay Assay to convert
#' @param reduction Name of DimReduc to set to main reducedDim in cds
#'
#' @rdname as.CellDataSet
#' @export
#' @method as.CellDataSet Seurat
#'
as.CellDataSet.Seurat <- function(x, assay = NULL, reduction = NULL, ...) {
  CheckDots(...)
  if (!PackageCheck('monocle', error = FALSE)) {
    stop("Please install monocle from Bioconductor before converting to a CellDataSet object")
  } else if (packageVersion(pkg = 'monocle') >= package_version(x = '2.99.0')) {
    stop("Seurat can only convert to/from Monocle v2.X objects")
  }
  assay <- assay %||% DefaultAssay(object = x)
  # make variables, then run `newCellDataSet`
  # create cellData counts
  counts <- GetAssayData(object = x, assay = assay, slot = "counts")
  # metadata
  cell.metadata <- x[[]]
  feature.metadata <- x[[assay]][[]]
  if (!"gene_short_name" %in% colnames(x = feature.metadata)) {
    feature.metadata$gene_short_name <- rownames(x = feature.metadata)
  }
  pd <- new(Class = "AnnotatedDataFrame", data = cell.metadata)
  fd <- new(Class = "AnnotatedDataFrame", data = feature.metadata)
  # Now, determine the expressionFamily
  if ("monocle" %in% names(x = Misc(object = x))) {
    expressionFamily <- Misc(object = x, slot = "monocle")[["expressionFamily"]]
  } else {
    if (all(counts == floor(x = counts))) {
      expressionFamily <- VGAM::negbinomial.size()
    } else if (any(counts < 0)) {
      expressionFamily <- VGAM::uninormal()
    } else {
      expressionFamily <- VGAM::tobit()
    }
  }
  cds <- monocle::newCellDataSet(
    cellData = counts,
    phenoData = pd,
    featureData = fd,
    expressionFamily = expressionFamily
  )
  if ("monocle" %in% names(x = Misc(object = x))) {
    monocle::cellPairwiseDistances(cds = cds) <- Misc(object = x, slot = "monocle")[["cellPairwiseDistances"]]
    monocle::minSpanningTree(cds = cds) <- Misc(object = x, slot = "monocle")[["minSpanningTree"]]
    Biobase::experimentData(cds = cds) <- Misc(object = x, slot = "monocle")[["experimentData"]]
    Biobase::protocolData(cds = cds) <- Misc(object = x, slot = "monocle")[["protocolData"]]
    Biobase::classVersion(cds = cds) <- Misc(object = x, slot = "monocle")[["classVersion"]]
    # no setter methods found for following slots
    slot(object = cds, name = "lowerDetectionLimit") <- Misc(object = x, slot = "monocle")[["lowerDetectionLimit"]]
    slot(object = cds, name = "dispFitInfo") <- Misc(object = x, slot = "monocle")[["dispFitInfo"]]
    slot(object = cds, name = "auxOrderingData") <- Misc(object = x, slot = "monocle")[["auxOrderingData"]]
    slot(object = cds, name = "auxClusteringData") <- Misc(object = x, slot = "monocle")[["auxClusteringData"]]
  }
  # adding dimensionality reduction data to the CDS
  dr.slots <- c("reducedDimS", "reducedDimK", "reducedDimW", "reducedDimA")
  reduction <- reduction %||% DefaultDimReduc(object = x, assay = assay)
  if (!is.null(x = reduction)) {
    if (grepl(pattern = 'tsne', x = tolower(x = reduction))) {
      slot(object = cds, name = "dim_reduce_type") <- "tSNE"
      monocle::reducedDimA(cds = cds) <- t(x = Embeddings(object = x[[reduction]]))
    } else {
      slot(object = cds, name = "dim_reduce_type") <- reduction
      monocle::reducedDimA(cds = cds) <- Loadings(object = x[[reduction]])
      slot(object = cds, name = "reducedDimS") <- Embeddings(object = x[[reduction]])
    }
    for (ii in dr.slots) {
      if (ii %in% names(x = slot(object = x[[reduction]], name = "misc"))) {
        slot(object = cds, name = ii) <- slot(object = x[[reduction]], name = "misc")[[ii]]
      }
    }
  }
  return(cds)
}

#' @rdname as.Graph
#' @export
#' @method as.Graph Matrix
#'
#' @examples
#' # converting sparse matrix
#' mat <- Matrix::rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
#' rownames(x = mat) <- paste0("feature_", 1:10)
#' colnames(x = mat) <- paste0("cell_", 1:10)
#' g <- as.Graph(x = mat)
#'
as.Graph.Matrix <- function(x, ...) {
  CheckDots(...)
  x <- as.sparse(x = x)
  if (is.null(x = rownames(x = x))) {
    stop("Please provide rownames to the matrix before converting to a Graph.")
  }
  if (is.null(x = colnames(x = x))) {
    stop("Please provide colnames to the matrix before converting to a Graph.")
  }
  return(as(object = x, Class = "Graph"))
}

#' @rdname as.Graph
#' @export
#' @method as.Graph matrix
#'
#' @examples
#' # converting dense matrix
#' mat <- matrix(data = 1:16, nrow = 4)
#' rownames(x = mat) <- paste0("feature_", 1:4)
#' colnames(x = mat) <- paste0("cell_", 1:4)
#' g <- as.Graph(x = mat)
#'
as.Graph.matrix <- function(x, ...) {
  CheckDots(...)
  return(as.Graph.Matrix(x = as(object = x, Class = 'Matrix')))
}

#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#'
#' @rdname as.sparse
#' @export
#' @method as.sparse data.frame
#'
as.sparse.data.frame <- function(x, ...) {
  CheckDots(...)
  return(as(object = as.matrix(x = x), Class = 'dgCMatrix'))
}

#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#'
#' @rdname as.sparse
#' @export
#' @method as.sparse Matrix
#'
as.sparse.Matrix <- function(x, ...) {
  CheckDots(...)
  return(as(object = x, Class = 'dgCMatrix'))
}

#' @rdname as.sparse
#' @export
#' @method as.sparse matrix
#'
as.sparse.matrix <- function(x, ...) {
  return(as.sparse.Matrix(x = x, ...))
}

#' @rdname Cells
#' @export
#'
Cells.default <- function(x) {
  return(colnames(x = x))
}

#' @rdname Cells
#' @export
#' @method Cells DimReduc
#'
Cells.DimReduc <- function(x) {
  return(rownames(x = x))
}

#' @param command Name of the command to pull, pass \code{NULL} to get the names of all commands run
#' @param value Name of the parameter to pull the value for
#'
#' @rdname Command
#' @export
#' @method Command Seurat
#'
Command.Seurat <- function(object, command = NULL, value = NULL, ...) {
  CheckDots(...)
  commands <- slot(object = object, name = "commands")
  if (is.null(x = command)) {
    return(names(x = commands))
  }
  if (is.null(x = commands[[command]])) {
    stop(command, " has not been run or is not a valid command.")
  }
  command <- commands[[command]]
  if (is.null(x = value)) {
    return(command)
  }
  params <- slot(object = command, name = "params")
  if (!value %in% names(x = params)) {
    stop(value, " is not a valid parameter for ", slot(object = command, name = "name"))
  }
  return(params[[value]])
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay DimReduc
#'
DefaultAssay.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'assay.used'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay Graph
#'
DefaultAssay.Graph <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.used'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay Seurat
#'
#' @examples
#' # Get current default assay
#' DefaultAssay(object = pbmc_small)
#'
DefaultAssay.Seurat <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'active.assay'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay SeuratCommand
#'
DefaultAssay.SeuratCommand <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.used'))
}

#' @export
#' @method DefaultAssay<- Assay
#'
"DefaultAssay<-.Assay" <- function(object, ..., value) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.used'))
  object <- UpdateSlots(object = object)
  slot(object = object, name = 'assay.orig') <- value
  return(object)
}

#' @export
#' @method DefaultAssay<- DimReduc
#'
"DefaultAssay<-.DimReduc" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'assay.used') <- value
  return(object)
}

#' @export
#' @method DefaultAssay<- Graph
#'
"DefaultAssay<-.Graph" <- function(object, ..., value) {
  object <- UpdateSlots(object = object)
  slot(object = object, name = 'assay.used') <- value
  return(object)
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay<- Seurat
#'
#' @examples
#' # Create dummy new assay to demo switching default assays
#' new.assay <- pbmc_small[["RNA"]]
#' Key(object = new.assay) <- "RNA2_"
#' pbmc_small[["RNA2"]] <- new.assay
#' # switch default assay to RNA2
#' DefaultAssay(object = pbmc_small) <- "RNA2"
#' DefaultAssay(object = pbmc_small)
#'
"DefaultAssay<-.Seurat" <- function(object, ..., value) {
  CheckDots(...)
  if (!value %in% names(x = slot(object = object, name = 'assays'))) {
    stop("Cannot find assay ", value)
  }
  slot(object = object, name = 'active.assay') <- value
  return(object)
}

#' @rdname Embeddings
#' @export
#' @method Embeddings DimReduc
#'
#' @examples
#' # Get the embeddings directly from a DimReduc object
#' Embeddings(object = pbmc_small[["pca"]])[1:5, 1:5]
#'
Embeddings.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'cell.embeddings'))
}

#' @param reduction Name of reduction to pull cell embeddings for
#'
#' @rdname Embeddings
#' @export
#' @method Embeddings Seurat
#'
#' @examples
#' # Get the embeddings from a specific DimReduc in a Seurat object
#' Embeddings(object = pbmc_small, reduction = "pca")[1:5, 1:5]
#'
Embeddings.Seurat <- function(object, reduction = 'pca', ...) {
  return(Embeddings(object = object[[reduction]], ...))
}

#' @param assay Assay to get
#'
#' @rdname GetAssay
#' @export
#' @method GetAssay Seurat
#'
#' @examples
#' GetAssay(object = pbmc_small, assay = "RNA")
#'
GetAssay.Seurat <- function(object, assay = NULL, ...) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  object.assays <- FilterObjects(object = object, classes.keep = 'Assay')
  if (!assay %in% object.assays) {
    stop(paste0(
      assay,
      " is not an assay present in the given object. Available assays are: ",
      paste(object.assays, collapse = ", ")
    ))
  }
  return(slot(object = object, name = 'assays')[[assay]])
}

#' @param slot Specific information to pull (i.e. counts, data, scale.data, ...)
#'
#' @rdname GetAssayData
#' @export
#' @method GetAssayData Assay
#'
#' @examples
#' # Get the data directly from an Assay object
#' GetAssayData(object = pbmc_small[["RNA"]], slot = "data")[1:5,1:5]
#'
GetAssayData.Assay <- function(object, slot = 'data', ...) {
  CheckDots(...)
  return(slot(object = object, name = slot))
}

#' @param assay Name of assay to pull data from
#'
#' @rdname GetAssayData
#' @export
#' @method GetAssayData Seurat
#'
#' @examples
#' # Get the data from a specific Assay in a Seurat object
#' GetAssayData(object = pbmc_small, assay = "RNA", slot = "data")[1:5,1:5]
#'
GetAssayData.Seurat <- function(object, slot = 'data', assay = NULL, ...) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  return(GetAssayData(
    object = GetAssay(object = object, assay = assay),
    slot = slot
  ))
}

#' @param assay Name of assay to pull highly variable feature information for
#'
#' @importFrom tools file_path_sans_ext
#'
#' @rdname HVFInfo
#' @export
#' @method HVFInfo Seurat
#'
#' @examples
#' # Get the HVF info from a specific Assay in a Seurat object
#' HVFInfo(object = pbmc_small, assay = "RNA")[1:5, ]
#'
HVFInfo.Seurat <- function(
  object,
  selection.method = NULL,
  assay = NULL,
  status = FALSE,
  ...
) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  if (is.null(x = selection.method)) {
    cmds <- apply(
      X = expand.grid(
        'FindVariableFeatures',
        FilterObjects(object = object, classes.keep = 'Assay')
      ),
      MARGIN = 1,
      FUN = paste,
      collapse = '.'
    )
    find.command <- Command(object = object)[Command(object = object) %in% cmds]
    if (length(x = find.command) < 1) {
      stop(
        "Please run 'FindVariableFeatures'",
        call. = FALSE
      )
    }
    find.command <- find.command[length(x = find.command)]
    test.command <- paste(file_path_sans_ext(x = find.command), assay, sep = '.')
    find.command <- ifelse(
      test = test.command %in% Command(object = object),
      yes = test.command,
      no = find.command
    )
    selection.method <- switch(
      EXPR = file_path_sans_ext(x = find.command),
      'FindVariableFeatures' = Command(
        object = object,
        command = find.command,
        value = 'selection.method'
      ),
      stop("Unknown command for finding variable features: '", find.command, "'", call. = FALSE)
    )
  }
  return(HVFInfo(
    object = GetAssay(object = object, assay = assay),
    selection.method = selection.method,
    status = status
  ))
}

#' @rdname Idents
#' @export
#' @method Idents Seurat
#'
Idents.Seurat <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'active.ident'))
}

#' @param cells Set cell identities for specific cells
#' @param drop Drop unused levels
#'
#' @rdname Idents
#' @export
#' @method Idents<- Seurat
#'
"Idents<-.Seurat" <- function(object, cells = NULL, drop = FALSE, ..., value) {
  CheckDots(...)
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  cells <- intersect(x = cells, y = colnames(x = object))
  cells <- match(x = cells, table = colnames(x = object))
  if (length(x = cells) == 0) {
    warning("Cannot find cells provided")
    return(object)
  }
  idents.new <- if (length(x = value) == 1 && value %in% colnames(x = object[[]])) {
    unlist(x = object[[value]], use.names = FALSE)[cells]
  } else {
    if (is.list(x = value)) {
      value <- unlist(x = value, use.names = FALSE)
    }
    rep_len(x = value, length.out = length(x = cells))
  }
  new.levels <- if (is.factor(x = idents.new)) {
    levels(x = idents.new)
  } else {
    unique(x = idents.new)
  }
  old.levels <- levels(x = object)
  levels <- c(new.levels, old.levels)
  idents.new <- as.vector(x = idents.new)
  idents <- as.vector(x = Idents(object = object))
  idents[cells] <- idents.new
  idents[is.na(x = idents)] <- 'NA'
  levels <- intersect(x = levels, y = unique(x = idents))
  names(x = idents) <- colnames(x = object)
  missing.cells <- which(x = is.na(x = names(x = idents)))
  if (length(x = missing.cells) > 0) {
    idents <- idents[-missing.cells]
  }
  idents <- factor(x = idents, levels = levels)
  slot(object = object, name = 'active.ident') <- idents
  if (drop) {
    object <- droplevels(x = object)
  }
  return(object)
}

#' @rdname IsGlobal
#' @export
#' @method IsGlobal default
#'
IsGlobal.default <- function(object, ...) {
  return(FALSE)
}

#' @rdname IsGlobal
#' @export
#' @method IsGlobal DimReduc
#'
IsGlobal.DimReduc <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'global'))
}

#' @rdname Key
#' @export
#' @method Key Assay
#'
#' @examples
#' # Get an Assay key
#' Key(object = pbmc_small[["RNA"]])
#'
Key.Assay <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'key'))
}

#' @rdname Key
#' @export
#' @method Key DimReduc
#'
#' @examples
#' # Get a DimReduc key
#' Key(object = pbmc_small[["pca"]])
#'
Key.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'key'))
}

#' @rdname Key
#' @export
#' @method Key Seurat
#'
#' @examples
#' # Show all keys associated with a Seurat object
#' Key(object = pbmc_small)
#'
Key.Seurat <- function(object, ...) {
  CheckDots(...)
  keyed.objects <- FilterObjects(object = object)
  return(sapply(
    X = keyed.objects,
    FUN = function(x) {
      return(Key(object = object[[x]]))
    }
  ))
}

#' @rdname Key
#' @export
#' @method Key<- Assay
#'
#' @examples
#' # Set the key for an Assay
#' Key(object = pbmc_small[["RNA"]]) <- "newkey_"
#' Key(object = pbmc_small[["RNA"]])
#'
"Key<-.Assay" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @rdname Key
#' @export
#' @method Key<- DimReduc
#'
#' @examples
#' # Set the key for DimReduc
#' Key(object = pbmc_small[["pca"]]) <- "newkey2_"
#' Key(object = pbmc_small[["pca"]])
#'
"Key<-.DimReduc" <- function(object, ..., value) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  old.key <- Key(object = object)
  slots <- Filter(
    f = function(x) {
      return(class(x = slot(object = object, name = x)) == 'matrix')
    },
    x = slotNames(x = object)
  )
  for (s in slots) {
    mat <- slot(object = object, name = s)
    if (!IsMatrixEmpty(x = mat)) {
      colnames(x = mat) <- sub(
        pattern = paste0('^', old.key),
        replacement = value,
        x = colnames(x = mat)
      )
    }
    slot(object = object, name = s) <- mat
  }
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @param projected Pull the projected feature loadings?
#'
#' @rdname Loadings
#' @export
#' @method Loadings DimReduc
#'
#' @examples
#' # Get the feature loadings for a given DimReduc
#' Loadings(object = pbmc_small[["pca"]])[1:5,1:5]
#'
Loadings.DimReduc <- function(object, projected = FALSE, ...) {
  CheckDots(...)
  projected <- projected %||% Projected(object = object)
  slot <- ifelse(
    test = projected,
    yes = 'feature.loadings.projected',
    no = 'feature.loadings'
  )
  return(slot(object = object, name = slot))
}

#' @param reduction Name of reduction to pull feature loadings for
#'
#' @rdname Loadings
#' @export
#' @method Loadings Seurat
#'
#' @examples
#' # Get the feature loadings for a specified DimReduc in a Seurat object
#' Loadings(object = pbmc_small, reduction = "pca")[1:5,1:5]
#'
Loadings.Seurat <- function(object, reduction = 'pca', projected = FALSE, ...) {
  return(Loadings(object = object[[reduction]], projected = projected, ...))
}

#' @rdname Loadings
#' @export
#' @method Loadings<- DimReduc
#'
#' @examples
#' # Set the feature loadings for a given DimReduc
#' new.loadings <- Loadings(object = pbmc_small[["pca"]])
#' new.loadings <- new.loadings + 0.01
#' Loadings(object = pbmc_small[["pca"]]) <- new.loadings
#'
"Loadings<-.DimReduc" <- function(object, projected = TRUE, ..., value) {
  CheckDots(...)
  slot.use <- ifelse(
    test = projected,
    yes = 'feature.loadings.projected',
    no = 'feature.loadings'
  )
  if (ncol(x = value) != length(x = object)) {
    stop("New feature loadings must have the dimensions as currently calculated")
  }
  slot(object = object, name = slot.use) <- value
  return(object)
}

#' @rdname Project
#' @export
#' @method Project Seurat
#'
Project.Seurat <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'project.name'))
}

#' @rdname Project
#' @export
#' @method Project<- Seurat
#'
"Project<-.Seurat" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'project.name') <- as.character(x = value)
  return(object)
}

#' @param slot Where to store the new data
#' @param new.data New data to insert
#'
#'
#' @importFrom stats na.omit
#'
#' @rdname SetAssayData
#' @export
#' @method SetAssayData Assay
#'
#' @examples
#' # Set an Assay slot directly
#' count.data <- GetAssayData(object = pbmc_small[["RNA"]], slot = "counts")
#' count.data <- as.matrix(x = count.data + 1)
#' new.assay <- SetAssayData(object = pbmc_small[["RNA"]], slot = "counts", new.data = count.data)
#'
SetAssayData.Assay <- function(object, slot, new.data, ...) {
  CheckDots(...)
  slots.use <- c('counts', 'data', 'scale.data')
  if (!slot %in% slots.use) {
    stop(
      "'slot' must be one of ",
      paste(slots.use, collapse = ', '),
      call. = FALSE
    )
  }
  if (!IsMatrixEmpty(x = new.data)) {
    if (any(grepl(pattern = '_', x = rownames(x = new.data)))) {
      warning(
        "Feature names cannot have underscores ('_'), replacing with dashes ('-')",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = new.data) <- gsub(
        pattern = '_',
        replacement = '-',
        x = rownames(x = new.data)
      )
    }
    if (ncol(x = new.data) != ncol(x = object)) {
      stop(
        "The new data doesn't have the same number of cells as the current data",
        call. = FALSE
      )
    }
    num.counts <- nrow(x = object)
    counts.names <- rownames(x = object)
    if (slot == 'scale.data' && nrow(x = new.data) > num.counts) {
      warning(
        "Adding more features than present in current data",
        call. = FALSE,
        immediate. = TRUE
      )
    } else if (slot %in% c('counts', 'data') && nrow(x = new.data) != num.counts) {
      warning(
        "The new data doesn't have the same number of features as the current data",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (!all(rownames(x = new.data) %in% counts.names)) {
      warning(
        "Adding features not currently present in the object",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    new.features <- na.omit(object = match(
      x = counts.names,
      table = rownames(x = new.data)
    ))
    new.cells <- colnames(x = new.data)
    if (!all(new.cells %in% colnames(x = object))) {
      stop(
        "All cell names must match current cell names",
        call. = FALSE
      )
    }
    new.data <- new.data[new.features, colnames(x = object), drop = FALSE]
    if (slot %in% c('counts', 'data') && !all(dim(x = new.data) == dim(x = object))) {
      stop(
        "Attempting to add a different number of cells and/or features",
        call. = FALSE
      )
    }
  }
  if (!is.vector(x = rownames(x = new.data))) {
    rownames(x = new.data) <- as.vector(x = rownames(x = new.data))
  }
  if (!is.vector(x = colnames(x = new.data))) {
    colnames(x = new.data) <- as.vector(x = colnames(x = new.data))
  }
  slot(object = object, name = slot) <- new.data
  return(object)
}

#' @param assay Name of assay whose data should be set
#'
#' @rdname SetAssayData
#' @export
#' @method SetAssayData Seurat
#'
#' @examples
#' # Set an Assay slot through the Seurat object
#' count.data <- GetAssayData(object = pbmc_small[["RNA"]], slot = "counts")
#' count.data <- as.matrix(x = count.data + 1)
#' new.seurat.object <- SetAssayData(
#'     object = pbmc_small,
#'     slot = "counts",
#'     new.data = count.data,
#'     assay = "RNA"
#' )
#'
SetAssayData.Seurat <- function(
  object,
  slot = 'data',
  new.data,
  assay = NULL,
  ...
) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  object[[assay]] <- SetAssayData(object = object[[assay]], slot = slot, new.data = new.data, ...)
  return(object)
}

#' @inheritParams Idents
#'
#' @rdname Idents
#' @export
#' @method SetIdent Seurat
#'
SetIdent.Seurat <- function(object, cells = NULL, value, ...) {
  #message(
  #  'With Seurat 3.X, setting identity classes can be done as follows:\n',
  #  'Idents(object = ',
  #  deparse(expr = substitute(expr = object)),
  #  if (!is.null(x = cells)) {
  #    paste0(', cells = ', deparse(expr = substitute(expr = cells)))
  #  },
  #  ') <- ',
  #  deparse(expr = substitute(expr = value))
  #)
  CheckDots(...)
  Idents(object = object, cells = cells) <- value
  return(object)
}

#' @inheritParams Idents
#' @param save.name Store current identity information under this name
#'
#' @rdname Idents
#' @export
#' @method StashIdent Seurat
#'
StashIdent.Seurat <- function(object, save.name = 'orig.ident', ...) {
  message(
    'With Seurat 3.X, stashing identity classes can be accomplished with the following:\n',
    deparse(expr = substitute(expr = object)),
    '[[',
    deparse(expr = substitute(expr = save.name)),
    ']] <- Idents(object = ',
    deparse(expr = substitute(expr = object)),
    ')'
  )
  CheckDots(...)
  object[[save.name]] <- Idents(object = object)
  return(object)
}

#' @rdname Stdev
#' @export
#' @method Stdev DimReduc
#'
#' @examples
#' # Get the standard deviations for each PC from the DimReduc object
#' Stdev(object = pbmc_small[["pca"]])
#'
Stdev.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'stdev'))
}

#' @param reduction Name of reduction to use
#'
#' @rdname Stdev
#' @export
#' @method Stdev Seurat
#'
#' @examples
#' # Get the standard deviations for each PC from the Seurat object
#' Stdev(object = pbmc_small, reduction = "pca")
#'
Stdev.Seurat <- function(object, reduction = 'pca', ...) {
  CheckDots(...)
  return(Stdev(object = object[[reduction]]))
}



#' @param assay Assay to subset on
#' @param ident.use Create a cell subset based on the provided identity classes
#' @param ident.remove Subtract out cells from these identity classes (used for
#' filtration)
#' @param max.cells.per.ident Can be used to downsample the data to a certain
#' max per cell ident. Default is INF.
#' @param random.seed Random seed for downsampling
#'
#' @rdname SubsetData
#' @export
#' @method SubsetData Seurat
#'
SubsetData.Seurat <- function(
  object,
  assay = NULL,
  cells = NULL,
  subset.name = NULL,
  ident.use = NULL,
  ident.remove = NULL,
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1,
  ...
) {
  .Deprecated(old = "SubsetData", new = "subset")
  expression <- character(length = 0L)
  if (!is.null(x = subset.name)) {
    sub <- gsub(
      pattern = '"',
      replacement = '',
      x = deparse(expr = substitute(expr = subset.name))
    )
    if (!is.infinite(x = low.threshold)) {
      expression <- c(
        expression,
        paste(sub, '>', deparse(expr = substitute(expr = low.threshold)))
      )
    }
    if (!is.infinite(x = high.threshold)) {
      expression <- c(
        expression,
        paste(sub, '<', deparse(expr = substitute(expr = high.threshold)))
      )
    }
    if (!is.null(x = accept.value)) {
      expression <- c(
        expression,
        paste(sub, '==', deparse(expr = substitute(expr = accept.value)))
      )
    }
  }
  #message(
  #  'With Seurat 3.X, subsetting Seurat objects can now be done with:\n',
  #  'subset(x = ',
  #  deparse(expr = substitute(expr = object)),
  #  if (length(x = expression) > 0) {
  #    paste0(', subset = ', paste(expression, collapse = ' & '))
  #  },
  #  if (length(x = c(cells, ident.use) > 0)) {
  #    paste0(', select = c("', paste0(c(cells, ident.use), collapse = '", '), '")')
  #  },
  #  if (!is.infinite(x = max.cells.per.ident)) {
  #    paste0(', downsample = ', max.cells.per.ident, ', seed = ', random.seed)
  #  },
  #  ')'
  #)
  assay <- assay %||% DefaultAssay(object = object)
  cells <- OldWhichCells(
    object = object,
    assay = assay,
    ident = ident.use,
    ident.remove = ident.remove,
    subset.name = subset.name,
    cells = cells,
    max.cells.per.ident = max.cells.per.ident,
    random.seed = random.seed,
    low.threshold = low.threshold,
    high.threshold = high.threshold,
    accept.value = accept.value,
    ...
  )
  # Subset all the Assays
  assays <- FilterObjects(object = object, classes.keep = 'Assay')
  for (assay in assays) {
    slot(object = object, name = "assays")[[assay]] <- SubsetData(
      object = object[[assay]],
      cells = cells
    )
  }
  # Subset all the DimReducs
  drs <- FilterObjects(object = object, classes.keep = 'DimReduc')
  for (dr in drs) {
    object[[dr]] <- CreateDimReducObject(
      embeddings = Embeddings(object = object[[dr]])[cells, ],
      loadings = Loadings(object = object[[dr]], projected = FALSE),
      projected = Loadings(object = object[[dr]], projected = TRUE),
      assay = DefaultAssay(object = object[[dr]]),
      stdev = Stdev(object = object[[dr]]),
      key = Key(object = object[[dr]]),
      jackstraw = slot(object = object[[dr]], name = "jackstraw"),
      misc = slot(object[[dr]], name = "misc")
    )
  }
  slot(object = object, name = "active.ident") <- Idents(object = object)[cells]
  slot(object = object, name = "meta.data") <- slot(object = object, name = "meta.data")[cells, ]
  return(object)
}

#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures Assay
#'
VariableFeatures.Assay <- function(object, selection.method = NULL, ...) {
  CheckDots(...)
  if (!is.null(x = selection.method)) {
    vf <- HVFInfo(object = object, selection.method = selection.method, status = TRUE)
    return(rownames(x = vf)[which(x = vf[, "variable"][, 1])])
  }
  return(slot(object = object, name = 'var.features'))
}

#' @param assay Name of assay to pull variable features for
#'
#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures Seurat
#'
VariableFeatures.Seurat <- function(object, assay = NULL, selection.method = NULL, ...) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  return(VariableFeatures(object = object[[assay]], selection.method = selection.method))
}

#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures<- Assay
#'
"VariableFeatures<-.Assay" <- function(object, ..., value) {
  CheckDots(...)
  if (length(x = value) == 0) {
    slot(object = object, name = 'var.features') <- character(length = 0)
    return(object)
  }
  if (any(grepl(pattern = '_', x = value))) {
    warning(
      "Feature names cannot have underscores '_', replacing with dashes '-'",
      call. = FALSE,
      immediate = TRUE
    )
    value <- gsub(pattern = '_', replacement = '-', x = value)
  }
  value <- split(x = value, f = value %in% rownames(x = object))
  if (length(x = value[['FALSE']]) > 0) {
    if (length(x = value[['TRUE']]) == 0) {
      stop("None of the features provided are in this Assay object", call. = FALSE)
    } else {
      warning(
        "Not all features provided are in this Assay object, removing the following feature(s): ",
        paste(value[['FALSE']], collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  slot(object = object, name = 'var.features') <- value[['TRUE']]
  return(object)
}

#' @inheritParams VariableFeatures.Seurat
#'
#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures<- Seurat
#'
"VariableFeatures<-.Seurat" <- function(object, assay = NULL, ..., value) {
  CheckDots(...)
  assay <- assay %||% DefaultAssay(object = object)
  VariableFeatures(object = object[[assay]]) <- value
  return(object)
}

#' @param idents A vector of identity classes to keep
#' @param slot Slot to pull feature data for
#' @param downsample Maximum number of cells per identity class, default is \code{Inf};
#' downsampling will happen after all other operations, including inverting the
#' cell selection
#' @param seed Random seed for downsampling. If NULL, does not set a seed
#'
#' @importFrom stats na.omit
#'
#' @rdname WhichCells
#' @export
#' @method WhichCells Seurat
#'
WhichCells.Seurat <- function(
  object,
  cells = NULL,
  idents = NULL,
  expression,
  slot = 'data',
  invert = FALSE,
  downsample = Inf,
  seed = 1,
  ...
) {
  CheckDots(...)
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  cell.order <- cells
  if (!is.null(x = idents)) {
    if (!is.null(x = seed)) {
      set.seed(seed = seed)
    }
    if (any(!idents %in% levels(x = Idents(object = object)))) {
      stop(
        "Cannot find the following identities in the object: ",
        paste(
          idents[!idents %in% levels(x = Idents(object = object))],
          sep = ', '
        )
      )
    }
    cells.idents <- unlist(x = lapply(
      X = idents,
      FUN = function(i) {
        cells.use <- which(x = as.vector(x = Idents(object = object)) == i)
        cells.use <- names(x = Idents(object = object)[cells.use])
        return(cells.use)
      }
    ))
    cells <- intersect(x = cells, y = cells.idents)
  }
  if (!missing(x = expression)) {
    objects.use <- FilterObjects(object = object)
    object.keys <- sapply(
      X = objects.use,
      FUN = function(i) {
        return(Key(object = object[[i]]))
      }
    )
    key.pattern <- paste0('^', object.keys, collapse = '|')
    expr <- if (is.call(x = substitute(expr = expression))) {
      substitute(expr = expression)
    } else {
      parse(text = expression)
    }
    expr.char <- as.character(x = expr)
    expr.char <- unlist(x = lapply(X = expr.char, FUN = strsplit, split = ' '))
    expr.char <- gsub(
      pattern = '(',
      replacement = '',
      x = expr.char,
      fixed = TRUE
    )
    expr.char <- gsub(
      pattern = '`',
      replacement = '',
      x = expr.char
    )
    vars.use <- which(
      x = expr.char %in% rownames(x = object) |
        expr.char %in% colnames(x = object[[]]) |
        grepl(pattern = key.pattern, x = expr.char, perl = TRUE)
    )
    data.subset <- FetchData(
      object = object,
      vars = expr.char[vars.use],
      cells = cells,
      slot = slot
    )
    data.subset <- subset.data.frame(x = data.subset, subset = eval(expr = expr))
    cells <- rownames(x = data.subset)
  }
  if (invert) {
    cell.order <- colnames(x = object)
    cells <- colnames(x = object)[!colnames(x = object) %in% cells]
  }
  cells <- CellsByIdentities(object = object, cells = cells)
  cells <- lapply(
    X = cells,
    FUN = function(x) {
      if (length(x = x) > downsample) {
        x <- sample(x = x, size = downsample, replace = FALSE)
      }
      return(x)
    }
  )
  cells <- as.character(x = na.omit(object = unlist(x = cells, use.names = FALSE)))
  cells <- cells[na.omit(object = match(x = cell.order, table = cells))]
  return(cells)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
#'
"$.Seurat" <- function(x, i, ...) {
  return(x[[i, drop = TRUE]])
}

#' @export
#'
"$.SeuratCommand" <- function(x, i, ...) {
  params <- slot(object = x, name = "params")
  return(params[[i]])
}

#' @export
#'
"$<-.Seurat" <- function(x, i, ..., value) {
  x[[i]] <- value
  return(x)
}

#' @export
#' @method [ Assay
#'
"[.Assay" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- 1:nrow(x = x)
  }
  if (missing(x = j)) {
    j <- 1:ncol(x = x)
  }
  return(GetAssayData(object = x)[i, j, ..., drop = FALSE])
}

#' @export
#' @method [ DimReduc
#'
"[.DimReduc" <- function(x, i, j, drop = FALSE, ...) {
  loadings <- Loadings(object = x)
  if (missing(x = i)) {
    i <- 1:nrow(x = loadings)
  }
  if (missing(x = j)) {
    j <- names(x = x)
  } else if (is.numeric(x = j)) {
    j <- names(x = x)[j]
  }
  bad.j <- j[!j %in% colnames(x = loadings)]
  j <- j[!j %in% bad.j]
  if (length(x = j) == 0) {
    stop("None of the requested loadings are present.")
  }
  if (length(x = bad.j) > 0) {
    warning(paste0("The following loadings are not present: ", paste(bad.j, collapse = ", ")))
  }
  return(Loadings(object = x)[i, j, drop = drop, ...])
}

#' @inheritParams subset.Seurat
#'
#' @rdname subset.Seurat
#' @export
#' @method [ Seurat
#'
#' @examples
#' pbmc_small[VariableFeatures(object = pbmc_small), ]
#' pbmc_small[, 1:10]
#'
"[.Seurat" <- function(x, i, j, ...) {
  if (missing(x = i) && missing(x = j)) {
    return(x)
  }
  if (missing(x = i)) {
    i <- NULL
  } else if (missing(x = j)) {
    j <- colnames(x = x)
  }
  if (is.logical(x = i)) {
    if (length(i) != nrow(x = x)) {
      stop("Incorrect number of logical values provided to subset features")
    }
    i <- rownames(x = x)[i]
  }
  if (is.logical(x = j)) {
    if (length(j) != ncol(x = x)) {
      stop("Incorrect number of logical values provided to subset cells")
    }
    j <- colnames(x = x)[j]
  }
  if (is.numeric(x = i)) {
    i <- rownames(x = x)[i]
  }
  if (is.numeric(x = j)) {
    j <- colnames(x = x)[j]
  }
  return(subset.Seurat(x = x, features = i, cells = j, ...))
}

#' @export
#' @method [ SeuratCommand
#'
"[.SeuratCommand" <- function(x, i, ...) {
  slot.use <- c("name", "timestamp", "call_string", "params")
  if (!i %in% slot.use) {
    stop("Invalid slot")
  }
  return(slot(object = x, name = i))
}

#' @export
#' @method [[ Assay
#'
"[[.Assay" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.features'))
  }
  data.return <- slot(object = x, name = 'meta.features')[, i, drop = FALSE, ...]
  if (drop) {
    data.return <- unlist(x = data.return, use.names = FALSE)
    names(x = data.return) <- rep.int(x = rownames(x = x), times = length(x = i))
  }
  return(data.return)
}

#' @export
#' @method [[ DimReduc
#'
"[[.DimReduc" <- function(x, i, j, drop = FALSE, ...) {
  if (missing(x = i)) {
    i <- 1:nrow(x = x)
  }
  if (missing(x = j)) {
    j <- names(x = x)
  } else if (is.numeric(x = j)) {
    j <- names(x = x)[j]
  }
  embeddings <- Embeddings(object = x)
  bad.j <- j[!j %in% colnames(x = embeddings)]
  j <- j[!j %in% bad.j]
  if (length(x = j) == 0) {
    stop("None of the requested embeddings are present.")
  }
  if (length(x = bad.j) > 0) {
    warning(paste0("The following embeddings are not present: ", paste(bad.j, collapse = ", ")))
  }
  return(embeddings[i, j, drop = drop, ...])
}

#' @export
#' @method [[ Seurat
#'
"[[.Seurat" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.data'))
  }
  if (length(x = i) == 0) {
    return(data.frame(row.names = colnames(x = x)))
  } else if (length(x = i) > 1 || any(i %in% colnames(x = slot(object = x, name = 'meta.data')))) {
    if (any(!i %in% colnames(x = slot(object = x, name = 'meta.data')))) {
      warning(
        "Cannot find the following bits of meta data: ",
        paste0(
          i[!i %in% colnames(x = slot(object = x, name = 'meta.data'))],
          collapse = ', '
        )
      )
    }
    i <- i[i %in% colnames(x = slot(object = x, name = 'meta.data'))]
    data.return <- slot(object = x, name = 'meta.data')[, i, drop = FALSE, ...]
    if (drop) {
      data.return <- unlist(x = data.return, use.names = FALSE)
      names(x = data.return) <- rep.int(x = colnames(x = x), times = length(x = i))
    }
  } else {
    slot.use <- unlist(x = lapply(
      X = c('assays', 'reductions', 'graphs', 'neighbors', 'commands'),
      FUN = function(s) {
        if (any(i %in% names(x = slot(object = x, name = s)))) {
          return(s)
        }
        return(NULL)
      }
    ))
    if (is.null(x = slot.use)) {
      stop("Cannot find '", i, "' in this Seurat object")
    }
    data.return <- slot(object = x, name = slot.use)[[i]]
  }
  return(data.return)
}

#' Coerce a SeuratCommand to a list
#'
#' @inheritParams base::as.list
#' @param complete Include slots besides just parameters (eg. call string, name, timestamp)
#'
#' @return A list with the parameters and, if \code{complete = TRUE}, the call string, name, and timestamp
#'
#' @export
#' @method as.list SeuratCommand
#'
as.list.SeuratCommand <- function(x, complete = FALSE, ...) {
  CheckDots(...)
  cmd <- slot(object = x, name = 'params')
  if (complete) {
    cmd <- append(
      x = cmd,
      values = sapply(
        X = grep(
          pattern = 'params',
          x = slotNames(x = x),
          invert = TRUE,
          value = TRUE
        ),
        FUN = slot,
        object = x,
        simplify = FALSE,
        USE.NAMES = TRUE
      ),
      after = 0
    )
  }
  for (i in 1:length(x = cmd)) {
    if (is.character(x = cmd[[i]])) {
      cmd[[i]] <- paste(trimws(x = cmd[[i]]), collapse = ' ')
    }
  }
  return(cmd)
}

#' @export
#' @method dim Assay
#'
dim.Assay <- function(x) {
  return(dim(x = GetAssayData(object = x)))
}

#' @export
#' @method dim DimReduc
#'
dim.DimReduc <- function(x) {
  return(dim(x = Embeddings(object = x)))
}

#' @export
#' @method dim Seurat
#'
dim.Seurat <- function(x) {
  return(dim(x = GetAssay(object = x)))
}

#' @export
#' @method dimnames Assay
#'
dimnames.Assay <- function(x) {
  return(dimnames(x = GetAssayData(object = x)))
}

#' @export
#' @method dimnames DimReduc
#'
dimnames.DimReduc <- function(x) {
  return(dimnames(x = Embeddings(object = x)))
}

#' @export
#' @method dimnames Seurat
#'
dimnames.Seurat <- function(x) {
  return(dimnames(x = GetAssay(object = x)))
}

#' @export
#' @method droplevels Seurat
#'
droplevels.Seurat <- function(x, ...) {
  slot(object = x, name = 'active.ident') <- droplevels(x = Idents(object = x), ...)
  return(x)
}

#' @export
#' @method length DimReduc
#'
length.DimReduc <- function(x) {
  return(ncol(x = Embeddings(object = x)))
}

#' @rdname Idents
#' @export
#' @method levels Seurat
#'
#' @examples
#' # Get the levels of identity classes of a Seurat object
#' levels(x = pbmc_small)
#'
levels.Seurat <- function(x) {
  return(levels(x = Idents(object = x)))
}

#' @rdname Idents
#' @export
#' @method levels<- Seurat
#'
#' @examples
#' # Reorder identity classes
#' levels(x = pbmc_small)
#' levels(x = pbmc_small) <- c('C', 'A', 'B')
#' levels(x = pbmc_small)
#'
"levels<-.Seurat" <- function(x, value) {
  idents <- Idents(object = x)
  if (!all(levels(x = idents) %in% value)) {
    stop("NA's generated by missing levels", call. = FALSE)
  }
  idents <- factor(x = idents, levels = value)
  Idents(object = x) <- idents
  return(x)
}



#' @export
#' @method names DimReduc
#'
names.DimReduc <- function(x) {
  return(colnames(x = Embeddings(object = x)))
}

#' @export
#' @method names Seurat
#'
names.Seurat <- function(x) {
  return(FilterObjects(object = x, classes.keep = c('Assay', 'DimReduc', 'Graph')))
}

#' Print the results of a dimensional reduction analysis
#'
#' Prints a set of features that most strongly define a set of components
#'
#' @param x An object
#' @param dims Number of dimensions to display
#' @param nfeatures Number of genes to display
#' @param projected Use projected slot
#' @param ... Arguments passed to other methods
#'
#' @return Set of features defining the components
#'
#' @aliases print
#' @seealso \code{\link[base]{cat}}
#'
#' @export
#' @method print DimReduc
#'
print.DimReduc <- function(x, dims = 1:5, nfeatures = 20, projected = FALSE, ...) {
  CheckDots(...)
  loadings <- Loadings(object = x, projected = projected)
  nfeatures <- min(nfeatures, nrow(x = loadings))
  if (ncol(x = loadings) == 0) {
    warning("Dimensions have not been projected. Setting projected = FALSE")
    projected <- FALSE
    loadings <- Loadings(object = x, projected = projected)
  }
  if (min(dims) > ncol(x = loadings)) {
    stop("Cannot print dimensions greater than computed")
  }
  if (max(dims) > ncol(x = loadings)) {
    warning(paste0("Only ", ncol(x = loadings), " dimensions have been computed."))
    dims <- min(dims):ncol(x = loadings)
  }
  for (dim in dims) {
    features <- TopFeatures(
      object = x,
      dim = dim,
      nfeatures = nfeatures * 2,
      projected = projected,
      balanced = TRUE
    )
    cat(Key(object = x), dim, '\n')
    pos.features <- split(x = features$positive, f = ceiling(x = seq_along(along.with = features$positive) / 10))
    cat("Positive: ", paste(pos.features[[1]], collapse = ", "), '\n')
    pos.features[[1]] <- NULL
    if (length(x = pos.features) > 0) {
      for (i in pos.features) {
        cat("\t  ", paste(i, collapse = ", "), '\n')
      }
    }
    neg.features <- split(x = features$negative, f = ceiling(x = seq_along(along.with = features$negative) / 10))
    cat("Negative: ", paste(neg.features[[1]], collapse = ", "), '\n')
    neg.features[[1]] <- NULL
    if (length(x = neg.features) > 0) {
      for (i in neg.features) {
        cat("\t  ", paste(i, collapse = ", "), '\n')
      }
    }
  }
}

#' Subset a Seurat object
#'
#' @param x Seurat object to be subsetted
#' @param subset Logical expression indicating features/variables to keep
#' @param i,features A vector of features to keep
#' @param j,cells A vector of cells to keep
#' @param idents A vector of identity classes to keep
#' @param ... Extra parameters passed to \code{\link{WhichCells}},
#' such as \code{slot}, \code{invert}, or \code{downsample}
#'
#' @return A subsetted Seurat object
#'
#' @rdname subset.Seurat
#' @aliases subset
#' @seealso \code{\link[base]{subset}} \code{\link{WhichCells}}
#'
#' @export
#' @method subset Seurat
#'
#' @examples
#' subset(x = pbmc_small, subset = MS4A1 > 4)
#' subset(x = pbmc_small, subset = `DLGAP1-AS1` > 2)
#' subset(x = pbmc_small, idents = '0', invert = TRUE)
#' subset(x = pbmc_small, subset = MS4A1 > 3, slot = 'counts')
#' subset(x = pbmc_small, features = VariableFeatures(object = pbmc_small))
#'
subset.Seurat <- function(x, subset, cells = NULL, features = NULL, idents = NULL, ...) {
  if (!missing(x = subset)) {
    subset <- deparse(expr = substitute(expr = subset))
  }
  cells <- WhichCells(
    object = x,
    cells = cells,
    idents = idents,
    expression = subset,
    ...
  )
  if (length(x = cells) == 0) {
    stop("No cells found", call. = FALSE)
  }
  if (all(cells %in% Cells(x = x)) && length(x = cells) == length(x = Cells(x = x)) && is.null(x = features)) {
    return(x)
  }
  assays <- FilterObjects(object = x, classes.keep = 'Assay')
  # Filter Assay objects
  for (assay in assays) {
    assay.features <- features %||% rownames(x = x[[assay]])
    slot(object = x, name = 'assays')[[assay]] <- tryCatch(
      expr = subset.Assay(x = x[[assay]], cells = cells, features = assay.features),
      error = function(e) {
        return(NULL)
      }
    )
  }
  slot(object = x, name = 'assays') <- Filter(
    f = Negate(f = is.null),
    x = slot(object = x, name = 'assays')
  )
  if (length(x = FilterObjects(object = x, classes.keep = 'Assay')) == 0 || is.null(x = x[[DefaultAssay(object = x)]])) {
    stop("Under current subsetting parameters, the default assay will be removed. Please adjust subsetting parameters or change default assay.", call. = FALSE)
  }
  # Filter DimReduc objects
  for (dimreduc in FilterObjects(object = x, classes.keep = 'DimReduc')) {
    x[[dimreduc]] <- tryCatch(
      expr = subset.DimReduc(x = x[[dimreduc]], cells = cells, features = features),
      error = function(e) {
        return(NULL)
      }
    )
  }
  # Remove metadata for cells not present
  slot(object = x, name = 'meta.data') <- slot(object = x, name = 'meta.data')[cells, , drop = FALSE]
  # Recalcualte nCount and nFeature
  for (assay in FilterObjects(object = x, classes.keep = 'Assay')) {
    n.calc <- CalcN(object = x[[assay]])
    if (!is.null(x = n.calc)) {
      names(x = n.calc) <- paste(names(x = n.calc), assay, sep = '_')
      x[[names(x = n.calc)]] <- n.calc
    }
  }
  # Filter metadata to keep nCount and nFeature for assays present
  ncolumns <- grep(
    pattern = '^nCount_|^nFeature_',
    x = colnames(x = x[[]]),
    value = TRUE
  )
  ncols.keep <- as.vector(x = outer(
    X = c('nCount_', 'nFeature_'),
    Y = FilterObjects(object = x, classes.keep = 'Assay'),
    FUN = paste0
  ))
  ncols.keep <- paste(ncols.keep, collapse = '|')
  ncols.remove <- grep(
    pattern = ncols.keep,
    x = ncolumns,
    value = TRUE,
    invert = TRUE
  )
  metadata.keep <- colnames(x = x[[]])[!colnames(x = x[[]]) %in% ncols.remove]
  slot(object = x, name = 'meta.data') <- x[[metadata.keep]]
  slot(object = x, name = 'graphs') <- list()
  Idents(object = x, drop = TRUE) <- Idents(object = x)[cells]
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname AddMetaData
#'
setMethod(
  f = '[[<-',
  signature = c('x' = 'Assay'),
  definition = function(x, i, ..., value) {
    meta.data <- x[[]]
    feature.names <- rownames(x = meta.data)
    if (is.data.frame(x = value)) {
      value <- lapply(
        X = 1:ncol(x = value),
        FUN = function(index) {
          v <- value[[index]]
          names(x = v) <- rownames(x = value)
          return(v)
        }
      )
    }
    err.msg <- "Cannot add more or fewer meta.features information without values being named with feature names"
    if (length(x = i) > 1) {
      # Add multiple bits of feature-level metadata
      value <- rep_len(x = value, length.out = length(x = i))
      for (index in 1:length(x = i)) {
        names.intersect <- intersect(x = names(x = value[[index]]), feature.names)
        if (length(x = names.intersect) > 0) {
          meta.data[names.intersect, i[index]] <- value[[index]][names.intersect]
        } else if (length(x = value) %in% c(nrow(x = meta.data), 1) %||% is.null(x = value)) {
          meta.data[i[index]] <- value[index]
        } else {
          stop(err.msg, call. = FALSE)
        }
      }
    } else {
      # Add a single column to feature-level metadata
      value <- unlist(x = value)
      if (length(x = intersect(x = names(x = value), y = feature.names)) > 0) {
        meta.data[, i] <- value[feature.names]
      } else if (length(x = value) %in% c(nrow(x = meta.data), 1) || is.null(x = value)) {
        meta.data[, i] <- value
      } else {
        stop(err.msg, call. = FALSE)
      }
    }
    slot(object = x, name = 'meta.features') <- meta.data
    return(x)
  }
)

#' @rdname AddMetaData
#'
setMethod( # because R doesn't allow S3-style [[<- for S4 classes
  f = '[[<-',
  signature = c('x' = 'Seurat'),
  definition = function(x, i, ..., value) {
    # Require names, no index setting
    if (!is.character(x = i)) {
      stop("'i' must be a character", call. = FALSE)
    }
    # Allow removing of other object
    if (is.null(x = value)) {
      slot.use <- if (i %in% colnames(x = x[[]])) {
        'meta.data'
      } else {
        FindObject(object = x, name = i)
      }
      if (is.null(x = slot.use)) {
        stop("Cannot find object ", i, call. = FALSE)
      }
      if (i == DefaultAssay(object = x)) {
        stop("Cannot delete the default assay", call. = FALSE)
      }
    }
    # remove disallowed characters from object name
    newi <- if (is.null(x = value)) {
      i
    } else {
      make.names(names = i)
    }
    if (any(i != newi)) {
      warning(
        "Invalid name supplied, making object name syntactically valid. New object name is ",
        newi,
        "; see ?make.names for more details on syntax validity",
        call. = FALSE,
        immediate. = TRUE
      )
      i <- newi
    }
    # Figure out where to store data
    slot.use <- if (inherits(x = value, what = 'Assay')) {
      # Ensure we have the same number of cells
      if (ncol(x = value) != ncol(x = x)) {
        stop(
          "Cannot add a different number of cells than already present",
          call. = FALSE
        )
      }
      # Ensure cell order stays the same
      if (all(Cells(x = value) %in% Cells(x = x)) && !all(Cells(x = value) == Cells(x = x))) {
        for (slot in c('counts', 'data', 'scale.data')) {
          assay.data <- GetAssayData(object = value, slot = slot)
          if (!IsMatrixEmpty(x = assay.data)) {
            assay.data <- assay.data[, Cells(x = x), drop = FALSE]
          }
          # Use slot because SetAssayData is being weird
          slot(object = value, name = slot) <- assay.data
        }
      }
      'assays'
    } else if (inherits(x = value, what = 'Graph')) {
      # Ensure Assay that Graph is associated with is present in the Seurat object
      if (is.null(x = DefaultAssay(object = value))) {
        warning(
          "Adding a Graph without an assay associated with it",
          call. = FALSE,
          immediate. = TRUE
        )
      } else if (!any(DefaultAssay(object = value) %in% Assays(object = x))) {
        stop("Cannot find assay '", DefaultAssay(object = value), "' in this Seurat object", call. = FALSE)
      }
      # Ensure Graph object is in order
      if (all(Cells(x = value) %in% Cells(x = x)) && !all(Cells(x = value) == Cells(x = x))) {
        value <- value[Cells(x = x), Cells(x = x)]
      }
      'graphs'
    } else if (inherits(x = value, what = 'DimReduc')) {
      # All DimReducs must be associated with an Assay
      if (is.null(x = DefaultAssay(object = value))) {
        stop("Cannot add a DimReduc without an assay associated with it", call. = FALSE)
      }
      # Ensure Assay that DimReduc is associated with is present in the Seurat object
      if (!IsGlobal(object = value) && !any(DefaultAssay(object = value) %in% Assays(object = x))) {
        stop("Cannot find assay '", DefaultAssay(object = value), "' in this Seurat object", call. = FALSE)
      }
      # Ensure DimReduc object is in order
      if (all(Cells(x = value) %in% Cells(x = x)) && !all(Cells(x = value) == Cells(x = x))) {
        slot(object = value, name = 'cell.embeddings') <- value[[Cells(x = x), ]]
      }
      'reductions'
    } else if (inherits(x = value, what = 'SeuratCommand')) {
      # Ensure Assay that SeuratCommand is associated with is present in the Seurat object
      if (is.null(x = DefaultAssay(object = value))) {
        warning(
          "Adding a command log without an assay associated with it",
          call. = FALSE,
          immediate. = TRUE
        )
      } else if (!any(DefaultAssay(object = value) %in% Assays(object = x))) {
        stop("Cannot find assay '", DefaultAssay(object = value), "' in this Seurat object", call. = FALSE)
      }
      'commands'
    } else if (is.null(x = value)) {
      slot.use
    } else {
      'meta.data'
    }
    if (slot.use == 'meta.data') {
      # Add data to object metadata
      meta.data <- x[[]]
      cell.names <- rownames(x = meta.data)
      # If we have metadata with names, ensure they match our order
      if (is.data.frame(x = value) && !is.null(x = rownames(x = value))) {
        meta.order <- match(x = rownames(x = meta.data), table = rownames(x = value))
        value <- value[meta.order, , drop = FALSE]
      }
      if (length(x = i) > 1) {
        # Add multiple pieces of metadata
        value <- rep_len(x = value, length.out = length(x = i))
        for (index in 1:length(x = i)) {
          meta.data[i[index]] <- value[index]
        }
      } else {
        # Add a single column to metadata
        if (length(x = intersect(x = names(x = value), y = cell.names)) > 0) {
          meta.data[, i] <- value[cell.names]
        } else if (length(x = value) %in% c(nrow(x = meta.data), 1) || is.null(x = value)) {
          meta.data[, i] <- value
        } else {
          stop("Cannot add more or fewer cell meta.data information without values being named with cell names", call. = FALSE)
        }
      }
      # Check to ensure that we aren't adding duplicate names
      if (any(colnames(x = meta.data) %in% FilterObjects(object = x))) {
        bad.cols <- colnames(x = meta.data)[which(colnames(x = meta.data) %in% FilterObjects(object = x))]
        stop(
          paste0(
            "Cannot add a metadata column with the same name as an Assay or DimReduc - ",
            paste(bad.cols, collapse = ", ")),
          call. = FALSE
        )
      }
      # Store the revised metadata
      slot(object = x, name = 'meta.data') <- meta.data
    } else {
      # Add other object to Seurat object
      # Ensure cells match in value and order
      if (!(class(x = value) %in% c('SeuratCommand', 'NULL')) && !all(Cells(x = value) == Cells(x = x))) {
        stop("All cells in the object being added must match the cells in this object", call. = FALSE)
      }
      # Ensure we're not duplicating object names
      if (!is.null(x = FindObject(object = x, name = i)) && !(class(x = value) %in% c(class(x = x[[i]]), 'NULL'))) {
        stop(
          "This object already contains ",
          i,
          " as a",
          ifelse(
            test = tolower(x = substring(text = class(x = x[[i]]), first = 1, last = 1)) %in% c('a', 'e', 'i', 'o', 'u'),
            yes = 'n ',
            no = ' '
          ),
          class(x = x[[i]]),
          "; duplicate names are not allowed",
          call. = FALSE
        )
      }
      # Check keyed objects
      if (inherits(x = value, what = c('Assay', 'DimReduc'))) {
        if (length(x = Key(object = value)) == 0) {
          Key(object = value) <- paste0(tolower(x = i), '_')
        }
        Key(object = value) <- UpdateKey(key = Key(object = value))
        # Check for duplicate keys
        object.keys <- sapply(
          X = FilterObjects(object = x),
          FUN = function(i) {
            return(Key(object = x[[i]]))
          }
        )
        if (Key(object = value) %in% object.keys && is.null(x = FindObject(object = x, name = i))) {
          # Attempt to create a duplicate key based off the name of the object being added
          new.keys <- c(paste0(tolower(x = i), c('_', paste0(RandomName(length = 2L), '_'))))
          # Select new key to use
          key.use <- min(which(x = !new.keys %in% object.keys))
          new.key <- if (is.infinite(x = key.use)) {
            RandomName(length = 17L)
          } else {
            new.keys[key.use]
          }
          warning(
            "Cannot add objects with duplicate keys (offending key: ",
            Key(object = value),
            "), setting key to '",
            new.key,
            "'",
            call. = FALSE
          )
          # Set new key
          Key(object = value) <- new.key
        }
      }
      # For Assays, run CalcN
      if (inherits(x = value, what = 'Assay')) {
        if ((!i %in% Assays(object = x)) |
            (i %in% Assays(object = x) && ! identical(
              x = GetAssayData(object = x, assay = i, slot = "counts"),
              y = GetAssayData(object = value, slot = "counts"))
            )) {
          n.calc <- CalcN(object = value)
          if (!is.null(x = n.calc)) {
            names(x = n.calc) <- paste(names(x = n.calc), i, sep = '_')
            x[[names(x = n.calc)]] <- n.calc
          }
        }
      }
      # When removing an Assay, clear out associated DimReducs, Graphs, and SeuratCommands
      if (is.null(x = value) && inherits(x = x[[i]], what = 'Assay')) {
        objs.assay <- FilterObjects(
          object = x,
          classes.keep = c('DimReduc', 'SeuratCommand', 'Graph')
        )
        objs.assay <- Filter(
          f = function(o) {
            return(all(DefaultAssay(object = x[[o]]) == i) && !IsGlobal(object = x[[o]]))
          },
          x = objs.assay
        )
        for (o in objs.assay) {
          x[[o]] <- NULL
        }
      }
      # If adding a command, ensure it gets put at the end of the command list
      if (inherits(x = value, what = 'SeuratCommand')) {
        slot(object = x, name = slot.use)[[i]] <- NULL
        slot(object = x, name = slot.use) <- Filter(
          f = Negate(f = is.null),
          x = slot(object = x, name = slot.use)
        )
      }
      slot(object = x, name = slot.use)[[i]] <- value
      slot(object = x, name = slot.use) <- Filter(
        f = Negate(f = is.null),
        x = slot(object = x, name = slot.use)
      )
    }
    CheckGC()
    return(x)
  }
)

setMethod(
  f = 'colMeans',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::colMeans(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colMeans',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(colMeans(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colSums',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::colSums(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colSums',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::colSums(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowMeans',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowMeans(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowMeans',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowMeans(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowSums',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowSums(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'rowSums',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowSums(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'show',
  signature = 'AnchorSet',
  definition = function(object) {
    cat('An AnchorSet object containing', nrow(x = slot(object = object, name = "anchors")),
        "anchors between", length(x = slot(object = object, name = "object.list")), "Seurat objects \n",
        "This can be used as input to IntegrateData or TransferData.")
  }
)

setMethod(
  f = 'show',
  signature = 'Assay',
  definition = function(object) {
    cat('Assay data with', nrow(x = object), 'features for', ncol(x = object), 'cells\n')
    if (length(x = VariableFeatures(object = object)) > 0) {
      top.ten <- head(x = VariableFeatures(object = object), n = 10L)
      top <- 'Top'
      variable <- 'variable'
    } else {
      top.ten <- head(x = rownames(x = object), n = 10L)
      top <- 'First'
      variable <- ''
    }
    features <- paste0(
      variable,
      ' feature',
      if (length(x = top.ten) != 1) {'s'}, ":\n"
    )
    features <- gsub(pattern = '^\\s+', replacement = '', x = features)
    cat(
      top,
      length(x = top.ten),
      features,
      paste(strwrap(x = paste(top.ten, collapse = ', ')), collapse = '\n'),
      '\n'
    )
  }
)

setMethod(
  f = 'show',
  signature = 'DimReduc',
  definition = function(object) {
    cat(
      "A dimensional reduction object with key", Key(object = object), '\n',
      'Number of dimensions:', length(x = object), '\n',
      'Projected dimensional reduction calculated: ', Projected(object = object), '\n',
      'Jackstraw run:', as.logical(x = JS(object = object)), '\n',
      'Computed using assay:', DefaultAssay(object = object), '\n'
    )
  }
)

setMethod(
  f = 'show',
  signature = 'JackStrawData',
  definition = function(object) {
    # empp <- GetJS(object = object, slot = "empirical.p.values")
    empp <- object$empirical.p.values
    # scored <- GetJS(object = object, slot = "overall.p.values")
    scored <- object$overall.p.values
    cat(
      "A JackStrawData object simulated on",
      nrow(x = empp),
      "features for",
      ncol(x = empp),
      "dimensions.\n",
      "Scored for:",
      nrow(x = scored),
      "dimensions.\n"
    )
  }
)

setMethod(
  f = "show",
  signature = "Seurat",
  definition = function(object) {
    assays <- FilterObjects(object = object, classes.keep = 'Assay')
    nfeatures <- sum(vapply(
      X = assays,
      FUN = function(x) {
        return(nrow(x = object[[x]]))
      },
      FUN.VALUE = integer(length = 1L)
    ))
    num.assays <- length(x = assays)
    cat("An object of class", class(x = object), "\n")
    cat(
      nfeatures,
      'features across',
      ncol(x = object),
      'samples within',
      num.assays,
      ifelse(test = num.assays == 1, yes = 'assay', no = 'assays'),
      "\n"
    )
    cat(
      "Active assay:",
      DefaultAssay(object = object),
      paste0('(', nrow(x = object), ' features)')
    )
    other.assays <- assays[assays != DefaultAssay(object = object)]
    if (length(x = other.assays) > 0) {
      cat(
        '\n',
        length(x = other.assays),
        'other',
        ifelse(test = length(x = other.assays) == 1, yes = 'assay', no = 'assays'),
        'present:',
        strwrap(x = paste(other.assays, collapse = ', '))
      )
    }
    reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
    if (length(x = reductions) > 0) {
      cat(
        '\n',
        length(x = reductions),
        'dimensional',
        ifelse(test = length(x = reductions) == 1, yes = 'reduction', no = 'reductions'),
        'calculated:',
        strwrap(x = paste(reductions, collapse = ', '))
      )
    }
    cat('\n')
  }
)

setMethod(
  f = 'show',
  signature = 'SeuratCommand',
  definition = function(object) {
    params <- slot(object = object, name = "params")
    params <- params[sapply(X = params, FUN = class) != "function"]
    cat(
      "Command: ", slot(object = object, name = "call.string"), '\n',
      "Time: ", as.character(slot(object = object, name = "time.stamp")), '\n',
      sep = ""
    )
    for(p in 1:length(params)){
      cat(
        names(params[p]), ":", params[[p]], "\n"
      )
    }
  }
)

setMethod(
  f = 'show',
  signature = 'seurat',
  definition = function(object) {
    cat(
      "An old seurat object\n",
      nrow(x = object@data),
      'genes across',
      ncol(x = object@data),
      'samples\n'
    )
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Internal AddMetaData defintion
#
# @param object An object
# @param metadata A vector, list, or data.frame with metadata to add
# @param col.name A name for meta data if not a named list or data.frame
#
# @return object with metadata added
#
.AddMetaData <- function(object, metadata, col.name = NULL) {
  if (is.null(x = col.name) && is.atomic(x = metadata)) {
    stop("'col.name' must be provided for atomic metadata types (eg. vectors)")
  }
  if (inherits(x = metadata, what = c('matrix', 'Matrix'))) {
    metadata <- as.data.frame(x = metadata)
  }
  col.name <- col.name %||% names(x = metadata) %||% colnames(x = metadata)
  if (is.null(x = col.name)) {
    stop("No metadata name provided and could not infer it from metadata object")
  }
  object[[col.name]] <- metadata
  # if (class(x = metadata) == "data.frame") {
  #   for (ii in 1:ncol(x = metadata)) {
  #     object[[colnames(x = metadata)[ii]]] <- metadata[, ii, drop = FALSE]
  #   }
  # } else {
  #   object[[col.name]] <- metadata
  # }
  return(object)
}

# Calculate nCount and nFeature
#
# @param object An Assay object
#
# @return A named list with nCount and nFeature
#
#' @importFrom Matrix colSums
#
CalcN <- function(object) {
  if (IsMatrixEmpty(x = GetAssayData(object = object, slot = "counts"))) {
    return(NULL)
  }
  return(list(
    nCount = colSums(x = object, slot = 'counts'),
    nFeature = colSums(x = GetAssayData(object = object, slot = 'counts') > 0)
  ))
}

# Get the names of objects within a Seurat object that are of a certain class
#
# @param object A Seurat object
# @param classes.keep A vector of names of classes to get
#
# @return A vector with the names of objects within the Seurat object that are of class \code{classes.keep}
#
FilterObjects <- function(object, classes.keep = c('Assay', 'DimReduc')) {
  slots <- na.omit(object = Filter(
    f = function(x) {
      sobj <- slot(object = object, name = x)
      return(is.list(x = sobj) && !is.data.frame(x = sobj) && !is.package_version(x = sobj))
    },
    x = slotNames(x = object)
  ))
  slots <- grep(pattern = 'tools', x = slots, value = TRUE, invert = TRUE)
  slots <- grep(pattern = 'misc', x = slots, value = TRUE, invert = TRUE)
  slots.objects <- unlist(
    x = lapply(
      X = slots,
      FUN = function(x) {
        return(names(x = slot(object = object, name = x)))
      }
    ),
    use.names = FALSE
  )
  object.classes <- sapply(
    X = slots.objects,
    FUN = function(i) {
      return(inherits(x = object[[i]], what = classes.keep))
    }
  )
  object.classes <- which(x = object.classes, useNames = TRUE)
  return(names(x = object.classes))
}

# Update a Key
#
# @param key A character to become a Seurat Key
#
# @return An updated Key that's valid for Seurat
#
UpdateKey <- function(key) {
  if (grepl(pattern = '^[[:alnum:]]+_', x = key)) {
    return(key)
  } else {
    new.key <-  gsub(
      pattern = "[[:^alnum:]]",
      replacement = "",
      x = key,
      perl = TRUE
    )
    new.key <- paste0(new.key, '_')
    if (new.key == '_') {
      new.key <- paste0(RandomName(length = 3), '_')
    }
    warning(
      "Keys should be one or more alphanumeric characters followed by an underscore, setting key from ",
      key,
      " to ",
      new.key,
      call. = FALSE,
      immediate. = TRUE
    )
    return(new.key)
  }
}

# Update slots in an object
#
# @param object An object to update
#
# @return \code{object} with the latest slot definitions
#
UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(mode = class(x = xobj), length = 1L)
    }
  }
  return(object)
}

# Find the collection of an object within a Seurat object
#
# @param object A Seurat object
# @param name Name of object to find
#
# @return The collection (slot) of the object
#
FindObject <- function(object, name) {
  collections <- c('assays', 'graphs', 'neighbors', 'reductions', 'commands')
  object.names <- lapply(
    X = collections,
    FUN = function(x) {
      return(names(x = slot(object = object, name = x)))
    }
  )
  names(x = object.names) <- collections
  object.names <- Filter(f = Negate(f = is.null), x = object.names)
  for (i in names(x = object.names)) {
    if (name %in% names(x = slot(object = object, name = i))) {
      return(i)
    }
  }
  return(NULL)
}

# Get the top
#
# @param data Data to pull the top from
# @param num Pull top \code{num}
# @param balanced Pull even amounts of from positive and negative values
#
# @return The top \code{num}
# @seealso \{code{\link{TopCells}}} \{code{\link{TopFeatures}}}
#
Top <- function(data, num, balanced) {
  top <- if (balanced) {
    num <- round(x = num / 2)
    data <- data[order(data, decreasing = TRUE), , drop = FALSE]
    positive <- head(x = rownames(x = data), n = num)
    negative <- rev(x = tail(x = rownames(x = data), n = num))
    list(positive = positive, negative = negative)
  } else {
    data <- data[rev(x = order(abs(x = data))), , drop = FALSE]
    top <- head(x = rownames(x = data), n = num)
    top[order(data[top, ])]
  }
  return(top)
}
