
#' Remove features based on different criteria
#'
#' Remove features with no variability (method = "constant"), missing values
#' (method = "missing"), low intensity (method = "zero") or identical peaks
#' (method = "identical.peaks").
#'
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object
#' @param assay Character or integer. Name or number of assay to be used for
#' filtering.
#' @param method Method to determine features to be removed: "constant",
#' "missing", "zero", "identical.peaks".
#' @param freq Numeric. If more or equal than freq*100 % of the samples fulfill
#' criterion the feature is removed (not used if method = "constant").
#' @param verbose Logical. Should number of removed features be reported?
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object with features removed
#'
#' @export

remove_features <- function(se,
                            assay,
                            method,
                            freq = 0.25,
                            verbose = FALSE) {

    if (is.character(assay) && !(assay %in% names(assays(se)))) {
        stop(paste("assay", assay, "not found!"))
    }
    if (freq < 0 | freq > 1) {
        stop("freq needs to be between 0 and 1 (inclusive)!")
    }

    expr = assays(se)[[assay]]

    if (method == "constant") {
        sd = apply(expr, 1, sd)
        ind.rm = which(sd == 0)
    } else if (method == "identical.peaks") {
        ind.rm = remove_features_identical_peaks(se = se,
                                                 cutoff = freq)
    } else {
        if (method == "missing") {
            crit = apply(expr, 1, function(x) {sum(is.na(x))}) / ncol(expr)
        } else if (method == "zero") {
            crit = apply(expr, 1, function(x) {sum(x == 0)}) / ncol(expr)
        } else {
            stop(paste("method", method, "not known!"))
        }

        ind.rm = which(crit >= freq)
    }
    if (length(ind.rm) > 0) {
        if (verbose) {
            print(paste(length(ind.rm), "features removed"))
        }
        se = se[-ind.rm, ]
    }
    return(se)

}


## partly based on code of function findCorrelation_fast in the R package caret
#' @keywords internal

remove_features_identical_peaks <- function(se,
                                            cutoff = 0.5) {

    mz = assays(se)$mz
    peak = assays(se)$peak.detection

    x = matrix(nrow = nrow(se),
               ncol = nrow(se))
    for (i in 1:(nrow(se) - 1)) {
        for (j in (i + 1):nrow(se)) {
            ind = which(peak[i, ] == 1 & peak[j, ] == 1)
            if (length(ind) > 0) {
                x[i, j] = sum(mz[i, ind] == mz[j, ind]) / length(ind)
            }
        }
    }

    average <- colMeans(abs(x), na.rm = TRUE)
    average <- as.numeric(as.factor(average))
    combsAboveCutoff <- which(abs(x) > cutoff)
    colsToCheck <- ceiling(combsAboveCutoff/nrow(x))
    rowsToCheck <- combsAboveCutoff%%nrow(x)
    colsToDiscard <- average[colsToCheck] > average[rowsToCheck]
    rowsToDiscard <- !colsToDiscard
    deletecol <- c(colsToCheck[colsToDiscard], rowsToCheck[rowsToDiscard])
    deletecol <- unique(deletecol)
    deletecol
}
