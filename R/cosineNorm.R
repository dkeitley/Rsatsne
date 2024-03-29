#' Cosine normalization
#'
#' Perform cosine normalization on the column vectors of an expression matrix.
#'
#' @param x A gene expression matrix with cells as columns and genes as rows.
#' @param mode A string specifying the output to be returned.
#' @param subset.row A vector specifying which features to use to compute the L2 norm.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization is to be performed.
#' Only used when \code{x} is a \linkS4class{DelayedArray} object.
#'
#' @details
#' Cosine normalization removes scaling differences between expression vectors.
#' In the context of batch correction, this is usually applied to remove differences between batches that are normalized separately.
#' For example, \code{\link{fastMNN}} uses this function on the log-expression vectors by default.
#'
#' Technically, separate normalization introduces scaling differences in the normalized expression, which should manifest as a shift in the log-transformed expression.
#' However, in practice, single-cell data will contain many small counts (where the log function is near-linear) or many zeroes (which remain zero when the pseudo-count is 1).
#' In these applications, scaling differences due to separate normalization are better represented as scaling differences in the log-transformed values.
#'
#' If applied to the raw count vectors, cosine normalization is similar to library size-related (i.e., L1) normalization.
#' However, we recommend using dedicated methods for computing size factors to normalize raw count data.
#'
#' While the default is to directly return the cosine-normalized matrix, it may occasionally be desirable to obtain the L2 norm,
#' e.g., to apply an equivalent normalization to other matrices.
#' This can be achieved by setting \code{mode} accordingly.
#'
#' The function will return a \linkS4class{DelayedMatrix} if \code{x} is a \linkS4class{DelayedMatrix}.
#' This aims to delay the calculation of cosine-normalized values for very large matrices.
#'
#' @return
#' If \code{mode="matrix"}, a double-precision matrix of the same dimensions as \code{X} is returned, containing cosine-normalized values.
#'
#' If \code{mode="l2norm"}, a double-precision vector is returned containing the L2 norm for each cell.
#'
#' If \code{mode="all"}, a named list is returned containing the fields \code{"matrix"} and \code{"l2norm"}, which are as described above.
#'
#' @author
#' Aaron Lun
#'
#' @seealso
#' \code{\link{mnnCorrect}} and \code{\link{fastMNN}}, where this function gets used.
#'
#' @examples
#' A <- matrix(rnorm(1000), nrow=10)
#' str(cosineNorm(A))
#' str(cosineNorm(A, mode="l2norm"))
#'
#' @export
#' @importFrom methods is
#' @importFrom Matrix colSums
#' @importFrom BiocParallel SerialParam
#' @importFrom DelayedArray setAutoBPPARAM getAutoBPPARAM
cosineNorm <- function(x, mode=c("matrix", "all", "l2norm"), subset.row=NULL, BPPARAM=SerialParam()) {
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    # Setting it up in the case of a parallelized colSums via DelayedArray.
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    l2 <- sqrt(colSums(x^2))

    mode <- match.arg(mode)
    if (mode=="l2norm") {
        return(l2)
    }

    mat <- .apply_cosine_norm(x, l2)
    if (mode=="matrix") {
        mat
    } else {
        list(matrix=mat, l2norm=l2)
    }
}

#' @importFrom scuttle normalizeCounts
.apply_cosine_norm <- function(x, l2) {
    l2 <- pmax(1e-8, l2) # protect against zero-L2.
    scuttle::normalizeCounts(x, size_factors=l2, center_size_factors=FALSE, log=FALSE)
}
