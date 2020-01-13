#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Normalize raw data
#'
#' Normalize count data per cell and transform to log scale
#'
#' @param data Matrix with the raw count data
#' @param scale.factor Scale the data. Default is 1e4
#' @param verbose Print progress
#'
#' @return Returns a matrix with the normalize and log transformed data
#'
#' @import Matrix
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' mat
#' mat_norm <- LogNormalize(data = mat)
#' mat_norm
#'
LogNormalize <- function(data, scale.factor = 1e4, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as(object = data, Class = "dgCMatrix")
  }
  # call Rcpp function to normalize
  if (verbose) {
    cat("Performing log-normalization\n", file = stderr())
  }
  norm.data <- LogNorm(data, scale_factor = scale.factor, display_progress = verbose)
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}

#' Normalize raw data to fractions
#'
#' Normalize count data to relative counts per cell by dividing by the total
#' per cell. Optionally use a scale factor, e.g. for counts per million (CPM)
#' use \code{scale.factor = 1e6}.
#'
#' @param data Matrix with the raw count data
#' @param scale.factor Scale the result. Default is 1
#' @param verbose Print progress
#' @return Returns a matrix with the relative counts
#'
#' @import Matrix
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' mat
#' mat_norm <- RelativeCounts(data = mat)
#' mat_norm
#'
RelativeCounts <- function(data, scale.factor = 1, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as(object = data, Class = "dgCMatrix")
  }
  if (verbose) {
    cat("Performing relative-counts-normalization\n", file = stderr())
  }
  norm.data <- data
  norm.data@x <- norm.data@x / rep.int(colSums(norm.data), diff(norm.data@p)) * scale.factor
  return(norm.data)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param selection.method How to choose top variable features. Choose one of :
#' \itemize{
#'   \item{vst:}{ First, fits a line to the relationship of log(variance) and
#'   log(mean) using local polynomial regression (loess). Then standardizes the
#'   feature values using the observed mean and expected variance (given by the
#'   fitted line). Feature variance is then calculated on the standardized values
#'   after clipping to a maximum (see clip.max parameter).}
#'   \item{mean.var.plot (mvp):}{ First, uses a function to calculate average
#'   expression (mean.function) and dispersion (dispersion.function) for each
#'   feature. Next, divides features into num.bin (deafult 20) bins based on
#'   their average expression, and calculates z-scores for dispersion within
#'   each bin. The purpose of this is to identify variable features while
#'   controlling for the strong relationship between variability and average
#'   expression.}
#'   \item{dispersion (disp):}{ selects the genes with the highest dispersion values}
#' }
#' @param loess.span (vst method) Loess span parameter used when fitting the
#' variance-mean relationship
#' @param clip.max (vst method) After standardization values larger than
#' clip.max will be set to clip.max; default is 'auto' which sets this value to
#' the square root of the number of cells
#' @param mean.function Function to compute x-axis value (average expression).
#'  Default is to take the mean of the detected (i.e. non-zero) values
#' @param dispersion.function Function to compute y-axis value (dispersion).
#' Default is to take the standard deviation of all values
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @param binning.method Specifies how the bins should be computed. Available
#' methods are:
#' \itemize{
#'   \item{equal_width:}{ each bin is of equal width along the x-axis [default]}
#'   \item{equal_frequency:}{ each bin contains an equal number of features (can
#'   increase statistical power to detect overdispersed features at high
#'   expression values, at the cost of reduced resolution along the x-axis)}
#' }
#' @param verbose show progress bar for calculations
#'
#' @rdname FindVariableFeatures
#' @export
#'
FindVariableFeatures.default <- function(
  object,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = 'auto',
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  if (!inherits(x = object, 'Matrix')) {
    object <- as(object = as.matrix(x = object), Class = 'Matrix')
  }
  if (!inherits(x = object, what = 'dgCMatrix')) {
    object <- as(object = object, Class = 'dgCMatrix')
  }
  if (selection.method == "vst") {
    if (clip.max == 'auto') {
      clip.max <- sqrt(x = ncol(x = object))
    }
    hvf.info <- data.frame(mean = rowMeans(x = object))
    hvf.info$variance <- SparseRowVar2(
      mat = object,
      mu = hvf.info$mean,
      display_progress = verbose
    )
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
    fit <- loess(
      formula = log10(x = variance) ~ log10(x = mean),
      data = hvf.info[not.const, ],
      span = loess.span
    )
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    # use c function to get variance after feature standardization
    hvf.info$variance.standardized <- SparseRowVarStd(
      mat = object,
      mu = hvf.info$mean,
      sd = sqrt(hvf.info$variance.expected),
      vmax = clip.max,
      display_progress = verbose
    )
    colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
  } else {
    if (!inherits(x = mean.function, what = 'function')) {
      stop("'mean.function' must be a function")
    }
    if (!inherits(x = dispersion.function, what = 'function')) {
      stop("'dispersion.function' must be a function")
    }
    feature.mean <- mean.function(object, verbose)
    feature.dispersion <- dispersion.function(object, verbose)
    names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = object)
    feature.dispersion[is.na(x = feature.dispersion)] <- 0
    feature.mean[is.na(x = feature.mean)] <- 0
    data.x.breaks <- switch(
      EXPR = binning.method,
      'equal_width' = num.bin,
      'equal_frequency' = c(
        -1,
        quantile(
          x = feature.mean[feature.mean > 0],
          probs = seq.int(from = 0, to = 1, length.out = num.bin)
        )
      ),
      stop("Unknown binning method: ", binning.method)
    )
    data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
    names(x = data.x.bin) <- names(x = feature.mean)
    mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
    sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = sd)
    feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
      sd.y[as.numeric(x = data.x.bin)]
    names(x = feature.dispersion.scaled) <- names(x = feature.mean)
    hvf.info <- data.frame(feature.mean, feature.dispersion, feature.dispersion.scaled)
    rownames(x = hvf.info) <- rownames(x = object)
    colnames(x = hvf.info) <- paste0('mvp.', c('mean', 'dispersion', 'dispersion.scaled'))
  }
  return(hvf.info)
}

#' @param nfeatures Number of features to select as top variable features;
#' only used when \code{selection.method} is set to \code{'dispersion'} or
#' \code{'vst'}
#' @param mean.cutoff A two-length numeric vector with low- and high-cutoffs for
#' feature means
#' @param dispersion.cutoff A two-length numeric vector with low- and high-cutoffs for
#' feature dispersions
#'
#' @rdname FindVariableFeatures
#' @export
#' @method FindVariableFeatures Assay
#'
FindVariableFeatures.Assay <- function(
  object,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = 'auto',
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 2000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE,
  ...
) {
  if (length(x = mean.cutoff) != 2 || length(x = dispersion.cutoff) != 2) {
    stop("Both 'mean.cutoff' and 'dispersion.cutoff' must be two numbers")
  }
  if (selection.method == "vst") {
    data <- GetAssayData(object = object, slot = "counts")
    # if (ncol(x = data) < 1 || nrow(x = data) < 1) {
    if (IsMatrixEmpty(x = data)) {
      warning("selection.method set to 'vst' but count slot is empty; will use data slot instead")
      data <- GetAssayData(object = object, slot = "data")
    }
  } else {
    data <- GetAssayData(object = object, slot = "data")
  }
  hvf.info <- FindVariableFeatures(
    object = data,
    selection.method = selection.method,
    loess.span = loess.span,
    clip.max = clip.max,
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    verbose = verbose,
    ...
  )
  object[[names(x = hvf.info)]] <- hvf.info
  hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 0), ]
  if (selection.method == "vst") {
    hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
  } else {
    hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE), , drop = FALSE]
  }
  selection.method <- switch(
    EXPR = selection.method,
    'mvp' = 'mean.var.plot',
    'disp' = 'dispersion',
    selection.method
  )
  top.features <- switch(
    EXPR = selection.method,
    'mean.var.plot' = {
      means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 1] < mean.cutoff[2])
      dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & (hvf.info[, 3] < dispersion.cutoff[2])
      rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    },
    'dispersion' = head(x = rownames(x = hvf.info), n = nfeatures),
    'vst' = head(x = rownames(x = hvf.info), n = nfeatures),
    stop("Unkown selection method: ", selection.method)
  )
  VariableFeatures(object = object) <- top.features
  vf.name <- ifelse(
    test = selection.method == 'vst',
    yes = 'vst',
    no = 'mvp'
  )
  vf.name <- paste0(vf.name, '.variable')
  object[[vf.name]] <- rownames(x = object[[]]) %in% top.features
  return(object)
}

#' @inheritParams FindVariableFeatures.Assay
#' @param assay Assay to use
#'
#' @rdname FindVariableFeatures
#' @export
#' @method FindVariableFeatures Seurat
#'
FindVariableFeatures.Seurat <- function(
  object,
  assay = NULL,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = 'auto',
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 2000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- FindVariableFeatures(
    object = assay.data,
    selection.method = selection.method,
    loess.span = loess.span,
    clip.max = clip.max,
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    nfeatures = nfeatures,
    mean.cutoff = mean.cutoff,
    dispersion.cutoff = dispersion.cutoff,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#'
#' @param normalization.method Method for normalization.
#'  \itemize{
#'   \item{LogNormalize: }{Feature counts for each cell are divided by the total
#'   counts for that cell and multiplied by the scale.factor. This is then
#'   natural-log transformed using log1p.}
#'   \item{CLR: }{Applies a centered log ratio transformation}
#'   \item{RC: }{Relative counts. Feature counts for each cell are divided by the total
#'   counts for that cell and multiplied by the scale.factor. No log-transformation is applied.
#'   For counts per million (CPM) set \code{scale.factor = 1e6}}
#' }
#' @param scale.factor Sets the scale factor for cell-level normalization
#' @param margin If performing CLR normalization, normalize across features (1) or cells (2)
# @param across If performing CLR normalization, normalize across either "features" or "cells".
#' @param block.size How many cells should be run in each chunk, will try to split evenly across threads
#' @param verbose display progress bar for normalization procedure
#'
#' @rdname NormalizeData
#' @export
#'
NormalizeData.default <- function(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  margin = 1,
  block.size = NULL,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  if (is.null(x = normalization.method)) {
    return(object)
  }
  normalized.data <- if (nbrOfWorkers() > 1) {
    norm.function <- switch(
      EXPR = normalization.method,
      'LogNormalize' = LogNormalize,
      'CLR' = CustomNormalize,
      'RC' = RelativeCounts,
      stop("Unknown normalization method: ", normalization.method)
    )
    if (normalization.method != 'CLR') {
      margin <- 2
    }
    tryCatch(
      expr = Parenting(parent.find = 'Seurat', margin = margin),
      error = function(e) {
        invisible(x = NULL)
      }
    )
    dsize <- switch(
      EXPR = margin,
      '1' = nrow(x = object),
      '2' = ncol(x = object),
      stop("'margin' must be 1 or 2")
    )
    chunk.points <- ChunkPoints(
      dsize = dsize,
      csize = block.size %||% ceiling(x = dsize / nbrOfWorkers())
    )
    normalized.data <- future_lapply(
      X = 1:ncol(x = chunk.points),
      FUN = function(i) {
        block <- chunk.points[, i]
        data <- if (margin == 1) {
          object[block[1]:block[2], , drop = FALSE]
        } else {
          object[, block[1]:block[2], drop = FALSE]
        }
        clr_function <- function(x) {
          return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
        }
        args <- list(
          data = data,
          scale.factor = scale.factor,
          verbose = FALSE,
          custom_function = clr_function, margin = margin
        )
        args <- args[names(x = formals(fun = norm.function))]
        return(do.call(
          what = norm.function,
          args = args
        ))
      }
    )
    do.call(
      what = switch(
        EXPR = margin,
        '1' = 'rbind',
        '2' = 'cbind',
        stop("'margin' must be 1 or 2")
      ),
      args = normalized.data
    )
  } else {
    switch(
      EXPR = normalization.method,
      'LogNormalize' = LogNormalize(
        data = object,
        scale.factor = scale.factor,
        verbose = verbose
      ),
      'CLR' = CustomNormalize(
        data = object,
        custom_function = function(x) {
          return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
        },
        margin = margin,
        verbose = verbose
        # across = across
      ),
      'RC' = RelativeCounts(
        data = object,
        scale.factor = scale.factor,
        verbose = verbose
      ),
      stop("Unkown normalization method: ", normalization.method)
    )
  }
  return(normalized.data)
}

#' @rdname NormalizeData
#' @export
#' @method NormalizeData Assay
#'
NormalizeData.Assay <- function(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  margin = 1,
  verbose = TRUE,
  ...
) {
  object <- SetAssayData(
    object = object,
    slot = 'data',
    new.data = NormalizeData(
      object = GetAssayData(object = object, slot = 'counts'),
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      verbose = verbose,
      margin = margin,
      ...
    )
  )
  return(object)
}

#' @param assay Name of assay to use
#'
#' @rdname NormalizeData
#' @export
#' @method NormalizeData Seurat
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' pmbc_small <- NormalizeData(object = pbmc_small)
#' }
#'
NormalizeData.Seurat <- function(
  object,
  assay = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  margin = 1,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- NormalizeData(
    object = assay.data,
    normalization.method = normalization.method,
    scale.factor = scale.factor,
    verbose = verbose,
    margin = margin,
    ...
  )
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}


#' @importFrom future nbrOfWorkers
#'
#' @param features Vector of features names to scale/center. Default is variable features.
#' @param vars.to.regress Variables to regress out (previously latent.vars in
#' RegressOut). For example, nUMI, or percent.mito.
#' @param latent.data Extra data to regress out, should be cells x latent data
#' @param split.by Name of variable in object metadata or a vector or factor defining
#' grouping of cells. See argument \code{f} in \code{\link[base]{split}} for more details
#' @param model.use Use a linear model or generalized linear model
#' (poisson, negative binomial) for the regression. Options are 'linear'
#' (default), 'poisson', and 'negbinom'
#' @param use.umi Regress on UMI count data. Default is FALSE for linear
#' modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to return for scaled data. The default is 10.
#' Setting this can help reduce the effects of feautres that are only expressed in
#' a very small number of cells. If regressing out latent variables and using a
#' non-linear model, the default is 50.
#' @param block.size Default size for number of feautres to scale at in a single
#' computation. Increasing block.size may speed up calculations but at an
#' additional memory cost.
#' @param min.cells.to.block If object contains fewer than this number of cells,
#' don't block for scaling calculations.
#' @param verbose Displays a progress bar for scaling procedure
#'
#' @importFrom future.apply future_lapply
#'
#' @rdname ScaleData
#' @export
#'
ScaleData.default <- function(
  object,
  features = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  features <- features %||% rownames(x = object)
  features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
  object <- object[features, , drop = FALSE]
  object.names <- dimnames(x = object)
  min.cells.to.block <- min(min.cells.to.block, ncol(x = object))
  suppressWarnings(expr = Parenting(
    parent.find = "ScaleData.Assay",
    features = features,
    min.cells.to.block = min.cells.to.block
  ))
  split.by <- split.by %||% TRUE
  split.cells <- split(x = colnames(x = object), f = split.by)
  CheckGC()
  if (!is.null(x = vars.to.regress)) {
    if (is.null(x = latent.data)) {
      latent.data <- data.frame(row.names = colnames(x = object))
    } else {
      latent.data <- latent.data[colnames(x = object), , drop = FALSE]
      rownames(x = latent.data) <- colnames(x = object)
    }
    if (any(vars.to.regress %in% rownames(x = object))) {
      latent.data <- cbind(
        latent.data,
        t(x = object[vars.to.regress[vars.to.regress %in% rownames(x = object)], ])
      )
    }
    # Currently, RegressOutMatrix will do nothing if latent.data = NULL
    if (verbose) {
      message("Regressing out ", paste(vars.to.regress, collapse = ', '))
    }
    chunk.points <- ChunkPoints(dsize = nrow(x = object), csize = block.size)
    if (nbrOfWorkers() > 1) { # TODO: lapply
      chunks <- expand.grid(
        names(x = split.cells),
        1:ncol(x = chunk.points),
        stringsAsFactors = FALSE
      )
      object <- future_lapply(
        X = 1:nrow(x = chunks),
        FUN = function(i) {
          row <- chunks[i, ]
          group <- row[[1]]
          index <- as.numeric(x = row[[2]])
          return(RegressOutMatrix(
            data.expr = object[chunk.points[1, index]:chunk.points[2, index], split.cells[[group]], drop = FALSE],
            latent.data = latent.data[split.cells[[group]], , drop = FALSE],
            features.regress = features,
            model.use = model.use,
            use.umi = use.umi,
            verbose = FALSE
          ))
        }
      )
      if (length(x = split.cells) > 1) {
        merge.indices <- lapply(
          X = 1:length(x = split.cells),
          FUN = seq.int,
          to = length(x = object),
          by = length(x = split.cells)
        )
        object <- lapply(
          X = merge.indices,
          FUN = function(x) {
            return(do.call(what = 'rbind', args = object[x]))
          }
        )
        object <- do.call(what = 'cbind', args = object)
      } else {
        object <- do.call(what = 'rbind', args = object)
      }
    } else {
      object <- lapply(
        X = names(x = split.cells),
        FUN = function(x) {
          if (verbose && length(x = split.cells) > 1) {
            message("Regressing out variables from split ", x)
          }
          return(RegressOutMatrix(
            data.expr = object[, split.cells[[x]], drop = FALSE],
            latent.data = latent.data[split.cells[[x]], , drop = FALSE],
            features.regress = features,
            model.use = model.use,
            use.umi = use.umi,
            verbose = verbose
          ))
        }
      )
      object <- do.call(what = 'cbind', args = object)
    }
    dimnames(x = object) <- object.names
    CheckGC()
  }
  if (verbose && (do.scale || do.center)) {
    msg <- paste(
      na.omit(object = c(
        ifelse(test = do.center, yes = 'centering', no = NA_character_),
        ifelse(test = do.scale, yes = 'scaling', no = NA_character_)
      )),
      collapse = ' and '
    )
    msg <- paste0(
      toupper(x = substr(x = msg, start = 1, stop = 1)),
      substr(x = msg, start = 2, stop = nchar(x = msg)),
      ' data matrix'
    )
    message(msg)
  }
  if (inherits(x = object, what = c('dgCMatrix', 'dgTMatrix'))) {
    scale.function <- FastSparseRowScale
  } else {
    object <- as.matrix(x = object)
    scale.function <- FastRowScale
  }
  if (nbrOfWorkers() > 1) {
    blocks <- ChunkPoints(dsize = length(x = features), csize = block.size)
    chunks <- expand.grid(
      names(x = split.cells),
      1:ncol(x = blocks),
      stringsAsFactors = FALSE
    )
    scaled.data <- future_lapply(
      X = 1:nrow(x = chunks),
      FUN = function(index) {
        row <- chunks[index, ]
        group <- row[[1]]
        block <- as.vector(x = blocks[, as.numeric(x = row[[2]])])
        data.scale <- scale.function(
          mat = object[features[block[1]:block[2]], split.cells[[group]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE
        )
        dimnames(x = data.scale) <- dimnames(x = object[features[block[1]:block[2]], split.cells[[group]]])
        suppressWarnings(expr = data.scale[is.na(x = data.scale)] <- 0)
        CheckGC()
        return(data.scale)
      }
    )
    if (length(x = split.cells) > 1) {
      merge.indices <- lapply(
        X = 1:length(x = split.cells),
        FUN = seq.int,
        to = length(x = scaled.data),
        by = length(x = split.cells)
      )
      scaled.data <- lapply(
        X = merge.indices,
        FUN = function(x) {
          return(suppressWarnings(expr = do.call(what = 'rbind', args = scaled.data[x])))
        }
      )
      scaled.data <- suppressWarnings(expr = do.call(what = 'cbind', args = scaled.data))
    } else {
      suppressWarnings(expr = scaled.data <- do.call(what = 'rbind', args = scaled.data))
    }
  } else {
    scaled.data <- matrix(
      data = NA_real_,
      nrow = nrow(x = object),
      ncol = ncol(x = object),
      dimnames = object.names
    )
    max.block <- ceiling(x = length(x = features) / block.size)
    for (x in names(x = split.cells)) {
      if (verbose) {
        if (length(x = split.cells) > 1 && (do.scale || do.center)) {
          message(gsub(pattern = 'matrix', replacement = 'from split ', x = msg), x)
        }
        pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
      }
      for (i in 1:max.block) {
        my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
        my.inds <- my.inds[my.inds <= length(x = features)]
        data.scale <- scale.function(
          mat = object[features[my.inds], split.cells[[x]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE
        )
        dimnames(x = data.scale) <- dimnames(x = object[features[my.inds], split.cells[[x]]])
        scaled.data[features[my.inds], split.cells[[x]]] <- data.scale
        rm(data.scale)
        CheckGC()
        if (verbose) {
          setTxtProgressBar(pb = pb, value = i)
        }
      }
      if (verbose) {
        close(con = pb)
      }
    }
  }
  dimnames(x = scaled.data) <- object.names
  scaled.data[is.na(x = scaled.data)] <- 0
  CheckGC()
  return(scaled.data)
}

#' @rdname ScaleData
#' @export
#' @method ScaleData Assay
#'
ScaleData.Assay <- function(
  object,
  features = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  use.umi <- ifelse(test = model.use != 'linear', yes = TRUE, no = use.umi)
  slot.use <- ifelse(test = use.umi, yes = 'counts', no = 'data')
  features <- features %||% VariableFeatures(object)
  if (length(x = features) == 0) {
    features <- rownames(x = GetAssayData(object = object, slot = slot.use))
  }
  object <- SetAssayData(
    object = object,
    slot = 'scale.data',
    new.data = ScaleData(
      object = GetAssayData(object = object, slot = slot.use),
      features = features,
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      split.by = split.by,
      model.use = model.use,
      use.umi = use.umi,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = scale.max,
      block.size = block.size,
      min.cells.to.block = min.cells.to.block,
      verbose = verbose,
      ...
    )
  )
  suppressWarnings(expr = Parenting(
    parent.find = "ScaleData.Seurat",
    features = features,
    min.cells.to.block = min.cells.to.block,
    use.umi = use.umi
  ))
  return(object)
}

#' @param assay Name of Assay to scale
#'
#' @rdname ScaleData
#' @export
#' @method ScaleData Seurat
#'
ScaleData.Seurat <- function(
  object,
  features = NULL,
  assay = NULL,
  vars.to.regress = NULL,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  if (any(vars.to.regress %in% colnames(x = object[[]]))) {
    latent.data <- object[[vars.to.regress[vars.to.regress %in% colnames(x = object[[]])]]]
  } else {
    latent.data <- NULL
  }
  if (is.character(x = split.by) && length(x = split.by) == 1) {
    split.by <- object[[split.by]]
  }
  assay.data <- ScaleData(
    object = assay.data,
    features = features,
    vars.to.regress = vars.to.regress,
    latent.data = latent.data,
    split.by = split.by,
    model.use = model.use,
    use.umi = use.umi,
    do.scale = do.scale,
    do.center = do.center,
    scale.max = scale.max,
    block.size = block.size,
    min.cells.to.block = min.cells.to.block,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}