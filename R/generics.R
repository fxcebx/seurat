#' Add in metadata associated with either cells or features.
#'
#' Adds additional data to the object. Can be any piece of information
#' associated with a cell (examples include read depth, alignment rate,
#' experimental batch, or subpopulation identity) or feature (ENSG name,
#' variance). To add cell level information, add to the Seurat object. If adding
#' feature-level metadata, add to the Assay object (e.g. object[["RNA"]]))
#'
#' @param x,object An object
#' @param i,col.name Name to store metadata or object as
#' @param value,metadata Metadata or object to add
#' @param j Ignored
#' @param ... Arguments passed to other methods
#'
#' @return An object with metadata or and object added
#'
#' @rdname AddMetaData
#' @export AddMetaData
#'
#' @aliases SeuratAccess
#'
#' @examples
#' cluster_letters <- LETTERS[Idents(object = pbmc_small)]
#' names(cluster_letters) <- colnames(x = pbmc_small)
#' pbmc_small <- AddMetaData(
#'   object = pbmc_small,
#'   metadata = cluster_letters,
#'   col.name = 'letter.idents'
#' )
#' head(x = pbmc_small[[]])
#'
AddMetaData <- function(object, metadata, col.name = NULL) {
  UseMethod(generic = 'AddMetaData', object = object)
}

#' Convert a matrix (or Matrix) to the Graph class.
#'
#' @param x The matrix to convert
#' @param ... Arguments passed to other methods (ignored for now)
#'
#' @rdname as.Graph
#' @export as.Graph
#'
as.Graph <- function(x, ...) {
  UseMethod(generic = "as.Graph", object = x)
}

#' Convert objects to Seurat objects
#'
#' @param x An object to convert to class \code{Seurat}
#' @param ... Arguments passed to other methods
#'
#' @rdname as.Seurat
#' @export as.Seurat
#'
as.Seurat <- function(x, ...) {
  UseMethod(generic = 'as.Seurat', object = x)
}

#' Convert between data frames and sparse matrices
#'
#' @param x An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{as.sparse}: A sparse representation of the input data
#'
#' @rdname as.sparse
#' @export as.sparse
#'
as.sparse <- function(x, ...) {
  UseMethod(generic = 'as.sparse', object = x)
}

#' Get cells present in an object
#'
#' @param x An object
#'
#' @return A vector of cell names
#'
#' @rdname Cells
#' @export Cells
#'
#' @examples
#' Cells(x = pbmc_small)
#'
Cells <- function(x) {
  UseMethod(generic = 'Cells', object = x)
}

#' Get SeuratCommands
#'
#' Pull information on previously run commands in the Seurat object.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Either a SeuratCommand object or the requested paramter value
#'
#' @rdname Command
#' @export Command
#'
Command <- function(object, ...) {
  UseMethod(generic = 'Command', object = object)
}

#' Get and set the default assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
  #' @return The name of the default assay
#'
#' @rdname DefaultAssay
#' @export DefaultAssay
#'
DefaultAssay <- function(object, ...) {
  UseMethod(generic = 'DefaultAssay', object = object)
}

#' @inheritParams DefaultAssay
#' @param value Name of assay to set as default
#'
#' @return An object with the new default assay
#'
#' @rdname DefaultAssay
#' @export DefaultAssay<-
#'
"DefaultAssay<-" <- function(object, ..., value) {
  UseMethod(generic = 'DefaultAssay<-', object = object)
}

#' Get cell embeddings
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Embeddings
#' @export Embeddings
#'
Embeddings <- function(object, ...) {
  UseMethod(generic = 'Embeddings', object = object)
}

#' Cluster Determination
#'
#' Identify clusters of cells by a shared nearest neighbor (SNN) modularity
#' optimization based clustering algorithm. First calculate k-nearest neighbors
#' and construct the SNN graph. Then optimize the modularity function to
#' determine clusters. For a full description of the algorithms, see Waltman and
#' van Eck (2013) \emph{The European Physical Journal B}. Thanks to Nigel
#' Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity
#' optimizer code in Rcpp!
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a Seurat object where the idents have been updated with new cluster info;
#' latest clustering results will be stored in object metadata under 'seurat_clusters'.
#' Note that 'seurat_clusters' will be overwritten everytime FindClusters is run
#'
#' @export
#'
#' @rdname FindClusters
#' @export FindClusters
#'
FindClusters <- function(object, ...) {
  UseMethod(generic = 'FindClusters', object = object)
}

#' Gene expression markers of identity classes
#'
#' Finds markers (differentially expressed genes) for identity classes
#'
#' @param object An object
#' @param ... Arguments passed to other methods and to specific DE methods

#' @return data.frame with a ranked list of putative markers as rows, and associated
#' statistics as columns (p-values, ROC score, etc., depending on the test used (\code{test.use})). The following columns are always present:
#' \itemize{
#'   \item \code{avg_logFC}: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
#'   \item \code{pct.1}: The percentage of cells where the gene is detected in the first group
#'   \item \code{pct.2}: The percentage of cells where the gene is detected in the second group
#'   \item \code{p_val_adj}: Adjusted p-value, based on bonferroni correction using all genes in the dataset
#' }
#'
#' @details p-value adjustment is performed using bonferroni correction based on
#' the total number of genes in the dataset. Other correction methods are not
#' recommended, as Seurat pre-filters genes using the arguments above, reducing
#' the number of tests performed. Lastly, as Aaron Lun has pointed out, p-values
#' should be interpreted cautiously, as the genes used for clustering are the
#' same genes tested for differential expression.
#'
#' @references McDavid A, Finak G, Chattopadyay PK, et al. Data exploration,
#' quality control and testing in single-cell qPCR-based gene expression experiments.
#' Bioinformatics. 2013;29(4):461-467. doi:10.1093/bioinformatics/bts714
#' @references Trapnell C, et al. The dynamics and regulators of cell fate
#' decisions are revealed by pseudotemporal ordering of single cells. Nature
#' Biotechnology volume 32, pages 381-386 (2014)
#' @references Andrew McDavid, Greg Finak and Masanao Yajima (2017). MAST: Model-based
#' Analysis of Single Cell Transcriptomics. R package version 1.2.1.
#' https://github.com/RGLab/MAST/
#' @references Love MI, Huber W and Anders S (2014). "Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2." Genome Biology.
#' https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#'
#' @export
#'
#' @examples
#' # Find markers for cluster 2
#' markers <- FindMarkers(object = pbmc_small, ident.1 = 2)
#' head(x = markers)
#'
#' # Take all cells in cluster 2, and find markers that separate cells in the 'g1' group (metadata
#' # variable 'group')
#' markers <- FindMarkers(pbmc_small, ident.1 = "g1", group.by = 'groups', subset.ident = "2")
#' head(x = markers)
#'
#' # Pass 'clustertree' or an object of class phylo to ident.1 and
#' # a node to ident.2 as a replacement for FindMarkersNode
#' pbmc_small <- BuildClusterTree(object = pbmc_small)
#' markers <- FindMarkers(object = pbmc_small, ident.1 = 'clustertree', ident.2 = 5)
#' head(x = markers)
#'
#' @rdname FindMarkers
#' @export FindMarkers
#'
#' @aliases FindMarkersNode
#'
FindMarkers <- function(object, ...) {
  UseMethod(generic = 'FindMarkers', object = object)
}

#' SNN Graph Construction
#'
#' Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset. We
#' first determine the k-nearest neighbors of each cell. We use this knn graph
#' to construct the SNN graph by calculating the neighborhood overlap
#' (Jaccard index) between every cell and its k.param nearest neighbors.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns the object with object@@snn filled
#'
#' @examples
#' pbmc_small
#' # Compute an SNN on the gene expression level
#' pbmc_small <- FindNeighbors(pbmc_small, features = VariableFeatures(object = pbmc_small))
#'
#' # More commonly, we build the SNN on a dimensionally reduced form of the data
#' # such as the first 10 principle components.
#'
#' pbmc_small <- FindNeighbors(pbmc_small, reduction = "pca", dims = 1:10)
#'
#' @rdname FindNeighbors
#' @export FindNeighbors
#'
FindNeighbors <- function(object, ...) {
  UseMethod(generic = 'FindNeighbors', object = object)
}

#' Find variable features
#'
#' Identifies features that are outliers on a 'mean variability plot'.
#'
#' For the mean.var.plot method:
#' Exact parameter settings may vary empirically from dataset to dataset, and
#' based on visual inspection of the plot. Setting the y.cutoff parameter to 2
#' identifies features that are more than two standard deviations away from the
#' average dispersion within a bin. The default X-axis function is the mean
#' expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
#' calculations are not performed in log-space, but the results are reported in
#' log-space - see relevant functions for exact details.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname FindVariableFeatures
#' @export FindVariableFeatures
#'
#' @aliases FindVariableGenes
#'
FindVariableFeatures <- function(object, ...) {
  UseMethod(generic = 'FindVariableFeatures', object = object)
}

#' Get an Assay object from a given Seurat object.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns an Assay object
#'
#' @rdname GetAssay
#' @export GetAssay
#'
GetAssay <- function(object, ...) {
  UseMethod(generic = 'GetAssay', object = object)
}

#' General accessor function for the Assay class
#'
#' This function can be used to pull information from any of the slots in the Assay class. For
#' example, pull one of the data matrices("counts", "data", or "scale.data").
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns info from requested slot
#'
#' @rdname GetAssayData
#' @export GetAssayData
#'
GetAssayData <- function(object, ...) {
  UseMethod(generic = 'GetAssayData', object = object)
}

#' Get highly variable feature information
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return A dataframe with feature means, dispersion, and scaled dispersion
#'
#' @rdname HVFInfo
#' @export HVFInfo
#'
HVFInfo <- function(object, ...) {
  UseMethod(generic = 'HVFInfo', object = object)
}

#' Get, set, and manipulate an object's identity classes
#'
#' @param x,object An object
#' @param ... Arguments passed to other methods; for \code{RenameIdents}: named
#' arguments as \code{old.ident = new.ident}; for \code{ReorderIdent}: arguments
#' passed on to \code{\link{FetchData}}
#'
#' @return \code{Idents}: The cell identies
#'
#' @rdname Idents
#' @export Idents
#'
#' @examples
#' # Get cell identity classes
#' Idents(object = pbmc_small)
#'
Idents <- function(object, ... ) {
  UseMethod(generic = 'Idents', object = object)
}

#' @inheritParams Idents
#' @param value The name of the identites to pull from object metadata or the identities themselves
#'
#' @return \code{Idents<-}: An object with the cell identites changed
#'
#' @rdname Idents
#' @export Idents<-
#'
#' @examples
#' # Set cell identity classes
#' # Can be used to set identities for specific cells to a new level
#' Idents(object = pbmc_small, cells = 1:4) <- 'a'
#' head(x = Idents(object = pbmc_small))
#'
#' # Can also set idents from a value in object metadata
#' colnames(x = pbmc_small[[]])
#' Idents(object = pbmc_small) <- 'RNA_snn_res.1'
#' levels(x = pbmc_small)
#'
"Idents<-" <- function(object, ..., value) {
  UseMethod(generic = 'Idents<-', object = object)
}

#' Is an object global/persistent?
#'
#' Typically, when removing \code{Assay} objects from an \code{Seurat} object,
#' all associated objects (eg. \code{DimReduc}, \code{Graph}, and \code{SeuratCommand} objects)
#' are removed as well. If an associated object is marked as global/persistent,
#' the associated object will remain even if its original assay was deleted
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{TRUE} if the object is global/persistent otherwise \code{FALSE}
#'
#' @rdname IsGlobal
#' @export IsGlobal
#'
#' @examples
#' IsGlobal(pbmc_small[['pca']])
#'
IsGlobal <- function(object, ...) {
  UseMethod(generic = 'IsGlobal', object = object)
}


#' Get a key
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Key
#' @export Key
#'
Key <- function(object, ...) {
  UseMethod(generic = 'Key', object = object)
}

#' Set a key
#'
#' @inheritParams Key
#' @param value Key value
#'
#' @rdname Key
#' @export Key<-
#'
"Key<-" <- function(object, ..., value) {
  UseMethod(generic = 'Key<-', object = object)
}

#' Get feature loadings
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Loadings
#' @export Loadings
#'
Loadings <- function(object, ...) {
  UseMethod(generic = 'Loadings', object = object)
}

#' Add feature loadings
#'
#' @inheritParams Loadings
#' @param value Feature loadings to add
#'
#' @rdname Loadings
#' @export Loadings<-
#'
"Loadings<-" <- function(object, ..., value) {
  UseMethod(generic = 'Loadings<-', object = object)
}


#' Normalize Data
#'
#' Normalize the count data present in a given assay.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns object after normalization
#'
#' @rdname NormalizeData
#' @export NormalizeData
#'
NormalizeData <- function(object, ...) {
  UseMethod(generic = 'NormalizeData', object = object)
}




#' Run Principal Component Analysis
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{PrintPCAParams}.
#'
#' @param object An object
#' @param ... Arguments passed to other methods and IRLBA
#'
#' @return Returns Seurat object with the PCA calculation stored in the reductions slot
#'
#' @export
#'
#' @rdname RunPCA
#' @export RunPCA
#'
RunPCA <- function(object, ...) {
  UseMethod(generic = 'RunPCA', object = object)
}

#' Run t-distributed Stochastic Neighbor Embedding
#'
#' Run t-SNE dimensionality reduction on selected features. Has the option of
#' running in a reduced dimensional space (i.e. spectral tSNE, recommended),
#' or running based on a set of genes. For details about stored TSNE calculation
#' parameters, see \code{PrintTSNEParams}.
#'
#' @param object Seurat object
#' @param ... Arguments passed to other methods and to t-SNE call (most commonly used is perplexity)
#'
#' @rdname RunTSNE
#' @export RunTSNE
#'
RunTSNE <- function(object, ...) {
  UseMethod(generic = 'RunTSNE', object = object)
}

#' Get and set project information
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Project information
#'
#' @rdname Project
#' @export Project
#'
Project <- function(object, ...) {
  UseMethod(generic = 'Project', object = object)
}

#' @param value Project information to set
#'
#' @return An object with project information added
#'
#' @rdname Project
#' @export Project<-
#'
"Project<-" <- function(object, ..., value) {
  UseMethod(generic = 'Project<-', object = object)
}

#' Compute Jackstraw scores significance.
#'
#' Significant PCs should show a p-value distribution that is
#' strongly skewed to the left compared to the null distribution.
#' The p-value for each PC is based on a proportion test comparing the number
#' of features with a p-value below a particular threshold (score.thresh), compared with the
#' proportion of features expected under a uniform distribution of p-values.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a Seurat object
#'
#' @author Omri Wurtzel
#' @seealso \code{\link{JackStrawPlot}}
#'
#' @rdname ScoreJackStraw
#' @export ScoreJackStraw
#'
ScoreJackStraw <- function(object, ...) {
  UseMethod(generic = 'ScoreJackStraw', object = object)
}

#' @return \code{SetIdent}: An object with new identity classes set
#'
#' @rdname Idents
#' @export SetIdent
#'
#' @examples
#' # Set cell identity classes using SetIdent
#' cells.use <- WhichCells(object = pbmc_small, idents = '1')
#' pbmc_small <- SetIdent(object = pbmc_small, cells = cells.use, value = 'B')
#'
SetIdent <- function(object, ...) {
  UseMethod(generic = 'SetIdent', object = object)
}

#' Return a subset of the Seurat object
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Takes either a list of cells to use as a subset, or a
#' parameter (for example, a gene), to subset on.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @rdname SubsetData
#' @export SubsetData
#'
#' @examples
#' \dontrun{
#' pbmc1 <- SubsetData(object = pbmc_small, cells = colnames(x = pbmc_small)[1:40])
#' pbmc1
#' }
#'
SubsetData <- function(object, ...) {
  UseMethod(generic = 'SubsetData', object = object)
}

#' Convert objects to CellDataSet objects
#'
#' @param x An object to convert to class \code{CellDataSet}
#' @param ... Arguments passed to other methods
#'
#' @rdname as.CellDataSet
#' @export as.CellDataSet
#'
as.CellDataSet <- function(x, ...) {
  UseMethod(generic = 'as.CellDataSet', object = x)
}

#' Run UMAP
#'
#' Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional
#' reduction technique. To run, you must first install the umap-learn python
#' package (e.g. via \code{pip install umap-learn}). Details on this package can be
#' found here: \url{https://github.com/lmcinnes/umap}. For a more in depth
#' discussion of the mathematics underlying UMAP, see the ArXiv paper here:
#' \url{https://arxiv.org/abs/1802.03426}.
#'
#' @param object An object
#' @param ... Arguments passed to other methods and UMAP
#'
#' @return Returns a Seurat object containing a UMAP representation
#'
#' @references McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and
#' Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' # Run UMAP map on first 5 PCs
#' pbmc_small <- RunUMAP(object = pbmc_small, dims = 1:5)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction = 'umap')
#' }
#'
#' @rdname RunUMAP
#' @export RunUMAP
#'
RunUMAP <- function(object, ...) {
  UseMethod(generic = 'RunUMAP', object = object)
}

#' Scale and center the data.
#'
#' Scales and centers features in the dataset. If variables are provided in vars.to.regress,
#' they are individually regressed against each feautre, and the resulting residuals are
#' then scaled and centered.
#'
#' ScaleData now incorporates the functionality of the function formerly known
#' as RegressOut (which regressed out given the effects of provided variables
#' and then scaled the residuals). To make use of the regression functionality,
#' simply pass the variables you want to remove to the vars.to.regress parameter.
#'
#' Setting center to TRUE will center the expression for each feautre by subtracting
#' the average expression for that feautre. Setting scale to TRUE will scale the
#' expression level for each feautre by dividing the centered feautre expression
#' levels by their standard deviations if center is TRUE and by their root mean
#' square otherwise.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname ScaleData
#' @export ScaleData
#'
ScaleData <- function(object, ...) {
  UseMethod(generic = 'ScaleData', object = object)
}



#' Setter for multimodal data
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return object with the assay data set
#'
#' @rdname SetAssayData
#' @export SetAssayData
#'
SetAssayData <- function(object, ...) {
  UseMethod(generic = 'SetAssayData', object = object)
}

#' @return \code{StashIdent}: An object with the identities stashed
#'
#' @rdname Idents
#' @export StashIdent
#'
#' @examples
#' head(x = pbmc_small[[]])
#' pbmc_small <- StashIdent(object = pbmc_small, save.name = 'idents')
#' head(x = pbmc_small[[]])
#'
StashIdent <- function(object, save.name, ...) {
  UseMethod(generic = 'StashIdent', object = object)
}

#' Get the standard deviations for an object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Stdev
#' @export Stdev
#'
Stdev <- function(object, ...) {
  UseMethod(generic = 'Stdev', object = object)
}

#' Get and set variable feature information
#'
#' @param object An object
#' @param selection.method Method used to set variable features
#' @param ... Arguments passed to other methods
#'
#' @rdname VariableFeatures
#' @export VariableFeatures
#'
VariableFeatures <- function(object, ...) {
  UseMethod(generic = 'VariableFeatures', object = object)
}

#' @inheritParams VariableFeatures
#' @param value A character vector of variable features
#'
#' @rdname VariableFeatures
#' @export VariableFeatures<-
#'
"VariableFeatures<-" <- function(object, ..., value) {
  UseMethod(generic = 'VariableFeatures<-', object = object)
}

#' Identify cells matching certain criteria
#'
#' Returns a list of cells that match a particular set of criteria such as
#' identity class, high/low values for particular PCs, ect..
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return A vector of cell names
#'
#' @seealso \code{\link{FetchData}}
#' @rdname WhichCells
#' @export WhichCells
#'
#' @examples
#' WhichCells(object = pbmc_small, idents = 2)
#' WhichCells(object = pbmc_small, expression = MS4A1 > 3)
#' levels(x = pbmc_small)
#' WhichCells(object = pbmc_small, idents = c(1, 2), invert = TRUE)
#'
WhichCells <- function(object, ...) {
  UseMethod(generic = 'WhichCells', object = object)
}
