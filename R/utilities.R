#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate the standard deviation of logged values
#'
#' Calculate SD of logged values in non-log space (return answer in log-space)
#'
#' @param x A vector of values
#'
#' @return Returns the standard deviation in log-space
#'
#' @importFrom stats sd
#'
#' @export
#'
#' @examples
#' ExpSD(x = c(1, 2, 3))
#'
ExpSD <- function(x) {
  return(log1p(x = sd(x = expm1(x = x))))
}

#' Calculate the variance of logged values
#'
#' Calculate variance of logged values in non-log space (return answer in
#' log-space)
#'
#' @param x A vector of values
#'
#' @return Returns the variance in log-space
#'
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#' ExpVar(x = c(1, 2, 3))
#'
ExpVar <- function(x) {
  return(log1p(x = var(x = expm1(x = x))))
}


#' Apply a ceiling and floor to all values in a matrix
#'
#' @param data Matrix or data frame
#' @param min all values below this min value will be replaced with min
#' @param max all values above this max value will be replaced with max
#' @return Returns matrix after performing these floor and ceil operations
#' @export
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 20, prob = 0.2 ), nrow = 5)
#' mat
#' MinMax(data = mat, min = 4, max = 5)
#'
MinMax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

#' Calculate the percentage of all counts that belong to a given set of features
#'
#' This function enables you to easily calculate the percentage of all the counts belonging to a
#' subset of the possible features for each cell. This is useful when trying to compute the percentage
#' of transcripts that map to mitochondrial genes for example. The calculation here is simply the
#' column sum of the matrix present in the counts slot for features belonging to the set divided by
#' the column sum for all features times 100.
#'
#' @param object A Seurat object
#' @param pattern A regex pattern to match features against
#' @param features A defined feature set. If features provided, will ignore the pattern matching
#' @param col.name Name in meta.data column to assign. If this is not null, returns a Seurat object
#' with the proportion of the feature set stored in metadata.
#' @param assay Assay to use
#'
#' @return Returns a vector with the proportion of the feature set or if md.name is set, returns a
#' Seurat object with the proportion of the feature set stored in metadata.
#' @importFrom Matrix colSums
#' @export
#'
#' @examples
#' # Calculate the proportion of transcripts mapping to mitochondrial genes
#' # NOTE: The pattern provided works for human gene names. You may need to adjust depending on your
#' # system of interest
#' pbmc_small[["percent.mt"]] <- PercentageFeatureSet(object = pbmc_small, pattern = "^MT-")
#'
PercentageFeatureSet <- function(
  object,
  pattern = NULL,
  features = NULL,
  col.name = NULL,
  assay = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = features) && !is.null(x = pattern)) {
    warning("Both pattern and features provided. Pattern is being ignored.")
  }
  features <- features %||% grep(pattern = pattern, x = rownames(x = object[[assay]]), value = TRUE)
  percent.featureset <- colSums(x = GetAssayData(object = object, assay = assay, slot = "counts")[features, , drop = FALSE])/
    object[[paste0("nCount_", assay)]] * 100
  if (!is.null(x = col.name)) {
    object <- AddMetaData(object = object, metadata = percent.featureset, col.name = col.name)
    return(object)
  }
  return(percent.featureset)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

# Set a default value if an object is NOT null
#
# @param lhs An object to set if it's NOT null
# @param rhs The value to provide if x is NOT null
#
# @return lhs if lhs is null, else rhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}


# Check the use of ...
#
# @param ... Arguments passed to a function that fall under ...
# @param fxns A list/vector of functions or function names
#
# @return ...
#
# @importFrom utils argsAnywhere getAnywhere
#' @importFrom utils isS3stdGeneric methods argsAnywhere isS3method
#
# @examples
#
CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(
    X = fxns,
    FUN = function(x) {
      x <- tryCatch(
        expr = if (isS3stdGeneric(f = x)) {
          as.character(x = methods(generic.function = x))
        } else {
          x
        },
        error = function(...) {
          return(x)
        }
      )
      x <- if (is.character(x = x)) {
        sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
      } else if (length(x = x) <= 1) {
        list(x)
      }
      return(sapply(
        X = x,
        FUN = function(f) {
          return(names(x = formals(fun = f)))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(X = fxn.args, FUN = is.null, FUN.VALUE = logical(length = 1L))
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found")
  } else if (any(fxn.null)) {
    warning(
      "The following functions passed could not be found: ",
      paste(names(x = which(x = fxn.null)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(
        X = Filter(f = nchar, x = dfxns),
        FUN = function(x) {
          if (isS3method(method = x)) {
            x <- unlist(x = strsplit(x = x, split = '\\.'))
            x <- x[length(x = x) - 1L]
          }
          return(x)
        },
        FUN.VALUE = character(length = 1L)
      )
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(unique(x = fx), collapse = ', ')
      )
      if (any(nchar(x = dfxns) < 1)) {
        message(
          "In addition, there is/are ",
          length(x = Filter(f = Negate(f = nchar), x = dfxns)),
          " other function(s) that accept(s) the dots"
        )
      }
    } else {
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      msg <- paste0(
        "The following arguments are not used: ",
        paste(unused, collapse = ', ')
      )
      switch(
        EXPR = getOption(x = "Seurat.checkdots"),
        "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
        "stop" = stop(msg),
        "silent" = NULL,
        stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent")
      )
      unused.hints <- sapply(X = unused, FUN = OldParamHints)
      names(x = unused.hints) <- unused
      unused.hints <- na.omit(object = unused.hints)
      if (length(x = unused.hints) > 0) {
        message(
          "Suggested parameter: ",
          paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
          "\n"
        )
      }
    }
  }
}

# Check a list of objects for duplicate cell names
#
# @param object.list List of Seurat objects
# @param verbose Print message about renaming
# @param stop Error out if any duplicate names exist
#
# @return Returns list of objects with duplicate cells renamed to be unique
#
CheckDuplicateCellNames <- function(object.list, verbose = TRUE, stop = FALSE) {
  cell.names <- unlist(x = lapply(X = object.list, FUN = colnames))
  if (any(duplicated(x = cell.names))) {
    if (stop) {
      stop("Duplicate cell names present across objects provided.")
    }
    if (verbose) {
      warning("Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.")
    }
    object.list <- lapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        return(RenameCells(
          object = object.list[[x]],
          new.names = paste0(Cells(x = object.list[[x]]), "_", x)
        ))
      }
    )
  }
  return(object.list)
}

# Call gc() to perform garbage collection
#
CheckGC <- function() {
  if (getOption(x = "Seurat.memsafe")) {
    gc(verbose = FALSE)
  }
}

# Extract delimiter information from a string.
#
# Parses a string (usually a cell name) and extracts fields based on a delimiter
#
# @param string String to parse.
# @param field Integer(s) indicating which field(s) to extract. Can be a vector multiple numbers.
# @param delim Delimiter to use, set to underscore by default.
#
# @return A new string, that parses out the requested fields, and (if multiple), rejoins them with the same delimiter
#
# @export
#
# @examples
# ExtractField(string = 'Hello World', field = 1, delim = '_')
#
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}


# Check if a matrix is empty
#
# Takes a matrix and asks if it's empty (either 0x0 or 1x1 with a value of NA)
#
# @param x A matrix
#
# @return Whether or not \code{x} is empty
#
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

# Melt a data frame
#
# @param x A data frame
#
# @return A molten data frame
#
Melt <- function(x) {
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}

# Parenting parameters from one environment to the next
#
# This function allows one to modifiy a parameter in a parent environement
# The primary use of this is to ensure logging functions store correct parameters
# if they've been modified by a child function or method
#
# @param parent.find Regex pattern of name of parent function call to modify.
# For example, this can be the class name for a method that was dispatched previously
# @param ... Parameter names and values to parent; both name and value must be supplied
# in the standard \code{name = value} format; any number of name/value pairs can be specified
#
# @return No return, modifies parent environment directly
#
# @examples
# Parenting(parent.find = 'Seurat', features = features[features > 7])
#
Parenting <- function(parent.find = 'Seurat', ...) {
  calls <- as.character(x = sys.calls())
  calls <- lapply(
    X = strsplit(x = calls, split = '(', fixed = TRUE),
    FUN = '[',
    1
  )
  parent.index <- grep(pattern = parent.find, x = calls)
  if (length(x = parent.index) != 1) {
    warning(
      "Cannot find a parent environment called ",
      parent.find,
      immediate. = TRUE,
      call. = FALSE
    )
  } else {
    to.parent <- list(...)
    if (length(x = to.parent) == 0) {
      warning("Nothing to parent", immediate. = TRUE, call. = FALSE)
    } else if (is.null(x = names(x = to.parent))) {
      stop("All input must be in a key = value pair")
    } else if (length(x = Filter(f = nchar, x = names(x = to.parent))) != length(x = to.parent)) {
      stop("All inputs must be named")
    } else {
      parent.environ <- sys.frame(which = parent.index)
      for (i in 1:length(x = to.parent)) {
        parent.environ[[names(x = to.parent)[i]]] <- to.parent[[i]]
      }
    }
  }
}

# Generate a random name
#
# Make a name from randomly sampled lowercase letters,
# pasted together with no spaces or other characters
#
# @param length How long should the name be
# @param ... Extra parameters passed to sample
#
# @return A character with nchar == length of randomly sampled letters
#
# @seealso \code{\link{sample}}
#
RandomName <- function(length = 5L, ...) {
  CheckDots(..., fxns = 'sample')
  return(paste(sample(x = letters, size = length, ...), collapse = ''))
}
