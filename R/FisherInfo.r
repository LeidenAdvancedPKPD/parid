# parid package: Calculation of parameter identifiability
#                indicators for nonlinear mixed effects models.
#
# Copyright (C) 2023 LAP&P Consultants BV
# Contact: info@lapp.nl
#          www.lapp.nl
#          Archimedesweg 31
#          2333 CM Leiden, The Netherlands
# 
# This file is part of the parid package.
# The parid package is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# The parid package is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the parid package.
# If not, see <https://www.gnu.org/licenses/>. 


#--------------------------------------------------------------------------------------------------------
#------------------------------------------ Internal functions ------------------------------------------
#--------------------------------------------------------------------------------------------------------


#------------------------------------------ checkSetSquareAndNames ------------------------------------------
#' Check whether the input is a square matrix and set its row or column names
#'
#' Checks whether the input \code{mat} is a square matrix (or \code{NULL}), and whether it has row and/or column names.
#' If it has both row and column names, checks whether they are the same.
#' If it has one of the two, sets the other equal to it.
#' Returns the adapted matrix and a list of errors.
#'
#' @param mat The object to check
#'
#' @return  a list of two elements.
#'   The first is the matrix \code{mat} with row or column names adapted as described above
#'   (if there were no errors).
#'   The second is an error message (as a string), or \code{NULL} if there were no errors.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
checkSetSquareAndNames <-function(mat) {
  if(is.null(mat)) return(list(matrix(0, nrow = 0, ncol = 0), NULL))
  if(!is.matrix(mat)) return(list(mat, " not a matrix"))
  if(nrow(mat) != ncol(mat)) return(list(mat, " not a square matrix"))
  if(nrow(mat) == 0) return(list(matrix(0, nrow = 0, ncol = 0), NULL))
  rwn <- row.names(mat)
  cln <- colnames(mat)
  if (is.null(cln) & is.null(rwn)) return(list(mat, " has no row or column names"))
  if (is.null(rwn)) rwn <- cln
  if (is.null(cln)) cln <- rwn
  if (!isTRUE(all.equal(rwn, cln))) return(list(mat, " has non-matching row and column names"))
  row.names(mat) <- rwn
  colnames(mat) <- cln
  return(list(mat, NULL))
}


#------------------------------------------ checkSymm ------------------------------------------
#' Check whether a matrix is symmetric
#'
#' Checks whether the matrix \code{mat} is symmetric.
#'
#' @param mat    The matrix to check.
#'
#' @return  A string with an error message if not symmetric,
#'   and \code{NULL} if symmetric.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
checkSymm <- function(mat) {
  if(!isTRUE(all.equal(t(mat), mat))) return("not symmetric")
  return(NULL)
}


#------------------------------------------ removeEmptyRowCols ------------------------------------------
#' Removes empty rows and columns from a symmetric matrix
#'
#' Removes rows and columns from the input matrices \code{mat} and \code{varmat} if they contain only
#' zeroes in \code{mat} and only \code{FALSE} in \code{varmat}.
#' The input matrics are assumed to be symmetric.
#' This can be used to reduce variance-covariance matrices.
#'
#' @param mat    The numeric matrix to check.
#' @param varmat The boolean matrix to check.
#'
#' @return  A list of two elements, namely the adapted versions of \code{mat} and \code{varmat}.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
removeEmptyRowCols <- function(mat, varmat) {
  # Assume symmetry so can check by row:
  row0 <- apply(mat, MARGIN = 1, FUN = function(row) all(row == 0))
  rowF <- apply(varmat, MARGIN = 1, FUN = function(row) all(!row))
  row  <- ifelse(row0 & rowF, 0, 1)  # 0 if row can be removed (& operates elementwise), 1 if not
  ind  <- ifelse(row %*% t(row) == 0, FALSE, TRUE)   # Matrix with FALSE for elements that can be removed and TRUE for those that should stay.
  nams <- colnames(mat)[row == 1]
  nams <- list(nams, nams)
  return(list(matrix(mat[ind], nrow = sum(row), dimnames = nams), matrix(varmat[ind], nrow = sum(row), dimnames = nams)))
}


#------------------------------------------ checkPosDef ------------------------------------------
#' Check whether a matrix is positive definite
#'
#' Checks whether the matrix \code{mat} is positive definite.
#'
#' @param mat    The matrix to check.
#' @param tol    The check for positive definiteness tests whether all eigenvalues have a real part larger than \code{tol}.
#'
#' @return  A string with an error message if not positive definite or not square,
#'   and \code{NULL} if positive definite. The empty matrix is by definition positive definite .
#'
#' @author Martijn van Noort
#' 
#' @noRd 
checkPosDef <- function(mat, tol = 1e-8) {
  if (nrow(mat) != ncol(mat)) return("not square")
  if (nrow(mat) == 0) return(NULL)
  evRe <- Re(eigen(mat + t(mat), symmetric = TRUE, only.values = TRUE)[["values"]])
  if(!all(evRe > tol)) return("not positive definite")
  return(NULL)
}


#------------------------------------------ checkSetCovMat ------------------------------------------
#' Check variance-covariance matrices and return them in a standard format
#'
#' Performs several checks on the variance-covariance matrix inputs of the methods that calculate the FIM,
#' and return them in a standard format.
#' Can also be used to standardize the input to \code{\link{normalizeFim}}.
#' For the inputs \code{omega} and \code{sigma} it is checked that they are matrices 
#' that are square, symmetric, positive definite and with matching row and column names.
#' Those names are set equal if necessary.
#' For the inputs \code{varomega} and \code{varsigma} it is checked that they are symmetric matrices,
#' or single values from which symmetric matrices are constructed.
#' Their row and column names are set to match those of \code{omega} and \code{sigma}.
#' 
#' Superfluous rows and columns are removed from these matrices.
#' That is, any row (and corresponding column) with only zero variances and covariances in \code{omega} and
#' only \code{FALSE} in \code{varomega} will be removed from both matrices, as it would not affect the outcome.
#' Likewise for \code{sigma} and \code{varsigma}.
#' The resulting matrices \code{omega} and \code{sigma} should not both have dimension zero.
#'
#' @param omega     The variance-covariance matrix of the individual effects, provided as a numeric matrix.
#' @param sigma     The variance-covariance matrix of the residual effects, provided as a numeric matrix.
#' @param varomega  Indicators which elements of \code{omega} should be considered as variable, provided
#'   as a boolean matrix or as a single boolean or \code{NULL}.
#'   In case of a single boolean, it is assumed that this single value applies to every matrix element
#'   (also the off-diagonal ones).
#'   The value \code{NULL} means all nonzero elements of \code{omega} are TRUE, and the others are \code{FALSE}.
#' @param varsigma  Indicators which elements of \code{sigma} should be considered as variable, provided
#'   as a boolean matrix or as a single boolean or \code{NULL}.
#'   In case of a single boolean, it is assumed that this single value applies to every matrix element
#'   (also the off-diagonal ones).
#'   The value \code{NULL} means all nonzero elements of \code{sigma} are TRUE, and the others are \code{FALSE}.
#' @param caller    String containing the signature of the calling function, for use in error messages.
#'   A typical value would be 'FisherInfo::calcFimFromMatrix'.
#'
#' @return  \code{NULL} if there were errors. In this case an error message will be printed. 
#'   Otherwise a list of four elements.
#'   The first two are the variance-covariance matrices, i.e. \code{omega} and \code{sigma}.
#'   The second pair are the indicator matrices \code{varomega} and \code{varsigma}.
#'   All four are adapted from input as described above.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
checkSetCovMat <- function(omega, sigma, varomega, varsigma, caller) {
  # Check and set row and column names of omega and sigma:
  out <- checkSetSquareAndNames(omega)
  if (!is.null(out[[2]])) {
    processErrors(paste0(caller, ": 'omega'", out[[2]], ".\nExiting.\n"))
    return(NULL)
  }
  omega <- out[[1]]
  out <- checkSetSquareAndNames(sigma)
  if (!is.null(out[[2]])) {
    processErrors(paste0(caller, ": 'sigma'", out[[2]], ".\nExiting.\n"))
    return(NULL)
  }
  sigma <- out[[1]]
  # Check that omega and sigma are symmetric:
  check <- checkSymm(omega)
  if (!is.null(check)) {
    processErrors(paste0(caller, ": 'omega' ", check, ".\nExiting.\n"))
    return(NULL)
  }
  check <- checkSymm(sigma)
  if (!is.null(check)) {
    processErrors(paste0(caller, ": 'sigma' ", check, ".\nExiting.\n"))
    return(NULL)
  }
  if (ncol(omega) == 0) {
    # Handle empty matrix separately:
    varomega <- matrix(FALSE, nrow = 0, ncol = 0)
  } else {
    # Set varomega to matrix:
    if (is.null(varomega)) varomega <- omega != 0
    if (!is.matrix(varomega)) {
      dm <- nrow(omega)
      varomega <- matrix(rep(varomega[[1]], dm*dm), nrow = dm, ncol = dm)
    }
    # Set row and column names of varomega:
    dimnames(varomega) <- dimnames(omega)
    # Check that varomega is symmetric:
    check <- checkSymm(varomega)
    if (!is.null(check)) {
      processErrors(paste0(caller, ": 'varomega' ", check, ".\nExiting.\n"))
      return(NULL)
    }
    # Remove irrelevant columns and rows from omega and varomega:
    out <- removeEmptyRowCols(omega, varomega)
    omega <- out[[1]]
    varomega <- out[[2]]
    # Check that omega is positive definite:
    check <- checkPosDef(omega)
    if (!is.null(check)) {
      processErrors(paste0(caller, ": 'omega' ", check, ".\nExiting.\n"))
      return(NULL)
    }
  }
  if (ncol(sigma) == 0) {
    # Handle empty matrix separately:
    varsigma <- matrix(FALSE, nrow = 0, ncol = 0)
  } else {
    # Set varsigma to matrix:
    if (is.null(varsigma)) varsigma <- sigma != 0
    if (!is.matrix(varsigma)) {
      dm <- nrow(sigma)
      varsigma <- matrix(rep(varsigma[[1]], dm*dm), nrow = dm, ncol = dm)
    }
    # Set row and column names of varsigma:
    dimnames(varsigma) <- dimnames(sigma)
    # Check that varsigma is symmetric:
    check <- checkSymm(varsigma)
    if (!is.null(check)) {
      processErrors(paste0(caller, ": 'varsigma' ", check, ".\nExiting.\n"))
      return(NULL)
    }
    # Remove irrelevant columns and rows from sigma and varsigma:
    out <- removeEmptyRowCols(sigma, varsigma)
    sigma <- out[[1]]
    varsigma <- out[[2]]
    # Check that sigma is positive definite:
    check <- checkPosDef(sigma)
    if (!is.null(check)) {
      processErrors(paste0(caller, ": 'sigma' ", check, ".\nExiting.\n"))
      return(NULL)
    }
  }
  if (nrow(omega) == 0 & nrow(sigma) == 0) {
    processErrors(paste0(caller, ": 'omega' and 'sigma' are both empty.\nExiting.\n"))
    return(NULL)
  }
  return(list(omega, sigma, varomega, varsigma))
}


#------------------------------------------ checkSetDerivCols ------------------------------------------
#' Check that the matrix \code{df} contains all required derivatives and return their column names
#'
#' The FIM calculation requires all first derivatives "dy_d<p>" and all second derivatives "d2y_dp_q" or "d2y_dq_p",
#' where p can be any element of \code{vartheta} or any rowname of \code{omega} or \code{sigma},
#' and q can be any rowname of \code{omega} or \code{sigma}.
#'
#' This function checks that the matrix \code{df} contains all these derivatives, and returns their names in \code{df}.
#' It also checks that the derivative names are returned in canonical order, i.e., corresponding to the order in
#' omega and sigma.
#'
#' @param df        Data frame of first and second order variations as generated by \code{\link{calcVariations}}
#'   or \code{\link{calcVariationsFim}}, that is, with columns 't', 'i', 'y', 'dy_d<v1>' and
#'   'd2y_d<v1>_<v2>', where variables v1 and v2 are replaced by names.
#' @param omega     The variance-covariance matrix of the individual effects, provided as a numeric matrix.
#' @param sigma     The variance-covariance matrix of the residual effects, provided as a numeric matrix.
#' @param vartheta  Vector of names of population parameters for which derivatives are to be calculated.
#'   Should be a subset of the names appearing in the variational matrix, or \code{NULL} for none.
#' @param caller    String containing the signature of the calling function, for use in error messages.
#'   A typical value would be 'FisherInfo::calcFimFromMatrix'.
#'
#' @return  \code{NULL} if there were errors. In this case an error message will be printed. 
#'   Otherwise a list of five elements, named "theta", "eta", "eps", "thetaEta" and "thetaEps" containing
#'   the column names in \code{df} of the derivatives of the indicated type.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
checkSetDerivCols <- function(df, omega, sigma, vartheta, caller) {
  # List all eta and eps in omega and sigma:
  vareta <- row.names(omega)
  vareps <- row.names(sigma)
  
  # Identify columns in df of different types of derivatives
  nam <- list("theta"    = if(length(vartheta) > 0) paste0("dy_d", vartheta) else NULL,
              "eta"      = if(length(vareta) > 0) paste0("dy_d", vareta) else NULL,
              "eps"      = if(length(vareps) > 0) paste0("dy_d", vareps) else NULL,
              "thetaEta" = if(length(vartheta) > 0 && length(vareta) > 0) {
                c(rbind(paste0("d2y_d", rep(vartheta, each = length(vareta)), "_", vareta),
                        paste0("d2y_d", rep(vareta, each = length(vartheta)), "_", vartheta)))
                } else NULL,
              "thetaEps" = if(length(vartheta) > 0 && length(vareps) > 0) {
                c(rbind(paste0("d2y_d", rep(vartheta, each = length(vareps)), "_", vareps),
                        paste0("d2y_d", rep(vareps, each = length(vartheta)), "_", vartheta)))
                } else NULL
  )  # Note that c(rbind()) interleaves the vectors. Hence the order of omega elements will be preserved.
  nam <- lapply(nam, function(nams) { nams[nams %in% names(df)] })
  # Check that the number of columns is correct:
  check <- mapply(nams = nam, expnr = list(length(vartheta), length(vareta), length(vareps), length(vartheta) * length(vareta),
                                           length(vartheta) * length(vareps)),
                  function(nams, expnr) length(nams) == expnr)
  if (any(!check)) {
    check <- names(nam)[!check]
    processErrors(paste0(caller, ": 'df' does not contain all derivatives for ",
               paste0(check, collapse = ", "), ".\nExiting.\n"))
    return(NULL)
  }
  # Check that "thetaEta" and "thetaEps" contain the derivatives in the same order as omega and sigma, respectively:
  check <- all(unlist(mapply(nams = nam[c("thetaEta", "thetaEps")], orignams = list(vareta, vareps), FUN = function(nams, orignams) {
    all(unlist(lapply(vartheta, function(th) {
      candidates <- grep(paste0("_d", th, "_|_", th, "$"), nams, value = TRUE)  # columns from nams that should correspond to orignams
      all(unlist(mapply(cand = candidates, orignam = orignams, function(cand, orignam) grepl(paste0("_d", orignam, "_|_", orignam, "$"), cand))))
    })))
  })))
  if (!check) {
    processErrors(paste0(caller, ": 'nam' does not contain 'thetaEta' and 'thetaEps' derivatives in correct order.\nExiting.\n"))
    return(NULL)
  }
  return(nam)
}


#------------------------------------------ createOmegaSigmaCases ------------------------------------------
#' Create data frame indicating which derivatives of random effects should be taken
#'
#' Creates a data frame with on each row a case description of a derivative of a random effect,
#' given a matrix indicating the variability of elements of a variance-covariance matrix.
#' 
#' @param varmat    A named square boolean matrix indicating which elements of a variance-covariance matrix
#'   should be considered as variable. Row names should equal column names.
#' @param isOmega   \code{TRUE} if the matrix \code{varmat} corresponds to individual effects,
#'   \code{FALSE} if residual effects.
#'
#' @return  Data frame with cases.
#' It contains numeric columns "row" and "col" indicating the row and column index in \code{varmat} of the case,
#'   a string column "nam" providing a description of the case, and a boolean column "omega" equal to \code{isOmega}.
#'   The column index is never higher than the row index.
#'   If row = col, then "nam" is set to the parameter name (i.e., the row name of \code{varmat}), and otherwise it 
#'   is set to "(<par1>,<par2>)", where <par1> and <par2> are the parameter names of row and column, respectively.
#'   The data frame may be empty if there were no cases to list.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
createOmegaSigmaCases <- function(varmat, isOmega) {
  cases <- data.frame(which(varmat & lower.tri(varmat, diag = TRUE), arr.ind = TRUE))
  cases[, "nam"] <- if (nrow(cases) == 0) rep("", 0) else {
    ifelse(cases[, "row"] == cases[, "col"], row.names(varmat)[cases[, "row"]],
           paste0("(", row.names(varmat)[cases[, "row"]], ",", row.names(varmat)[cases[, "col"]], ")"))
  }
  cases[, "omega"] <- if (nrow(cases) == 0) rep(FALSE, 0) else isOmega
  return(cases)
}


#------------------------------------------ traceprod ------------------------------------------
#' Calculate the trace of the product of two given matrices
#'
#' The trace of the "skew" product of equal dimension matrices A and B is given by
#' tr(A^t B) = sum_ij A_ij B_ij.
#' 
#' @param mat1      A numeric matrix.
#' @param mat2      Another numeric matrix, of the same dimension.
#'
#' @return  tr(mat1^t mat2).
#'
#' @author Martijn van Noort
#' 
#' @noRd 
traceprod <- function(mat1, mat2) {
  # compute trace(mat1 * mat2) = sum_ij mat1[i, j]*mat2[i, j]
  sum(mat1 * mat2)
}


#------------------------------------------ computeNormalizers ------------------------------------------
#' Normalizing factors for variance matrix
#'
#' Computes normalizing factors for a given variance-covariance matrix.
#' The factor takes into account prior normalization.
#' This function is for internal use and performs no validation.
#'
#' @param mat        Variance-covariance matrix, as numeric matrix.
#'   Should be symmetric and positive definite.
#'   Columns and rows should be named, and the names should be the same.
#' @param relvec     Boolean vector of the same dimension as \code{mat},
#'   indicating which elements need to be normalized (\code{TRUE}) or not (\code{FALSE}).
#'   Elements need not be named.
#' @param existNorm  Named boolean vector specifying the existing normalization.
#'   This should have the same dimension as \code{mat}, and the names should equal the
#'   column names of \code{mat}, in the same order.
#'
#' @return   A data frame with variables 'row', 'col' and 'nam' specifying the matrix element,
#'   and 'value' containing the normalizing factor.
#'   This data frame may be empty if \code{mat} is empty.
#'   The function displays an error and returns \code{NULL} if \code{existNorm} is not in the specified format.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
computeNormalizers <- function(mat, relvec, existNorm) {
  if ((length(existNorm) != nrow(mat)) || any(names(existNorm) != colnames(mat))) {
    processErrors("FisherInfo::computeNormalizers: existing normalization information does not match variance-covariance matrix. Exiting.\n")
    return(NULL)
  }
  if (nrow(mat) == 0) {
    return(data.frame(row = character(0), col = character(0), value = numeric(0), nam = character(0)))
  }
  diagElts <- ifelse(relvec, sqrt(diag(mat)), 1) # Elements sqrt(mat) or 1 (if no normalization)
  fac <- ifelse(existNorm, sqrt(diag(mat)), 1)
  fac[fac <= 0] <- 1  # if fac == 0 then prior normalization would have normalized to 0, in which case the value of 'fac' does not matter.
  diagElts <- diagElts / fac # Undo prior normalization
  factors <- diagElts %*% t(diagElts)            # Matrix with multiplication factors for each element; factor = 1 if element should not be normalized.
  row.names(factors) <- colnames(factors) <- row.names(mat)
  cases <- reshape2::melt(cbind(data.frame(row = row.names(factors)), as.data.frame(factors)), id.vars = "row", variable.name = "col")
  cases[, "nam"] <- ifelse(cases[, "row"] == cases[, "col"], as.character(cases[, "row"]),
                           paste0("(", cases[, "row"], ",", cases[, "col"], ")"))
  return(cases)
}


#------------------------------------------ doCalcFim ------------------------------------------
#' Calculate Fisher Information Matrix
#'
#' This is the workhorse for all FIM calculation methods.
#' The difference with \code{calcFimFromMatrix} is that it does not check the validity of the 
#' variance-covariance matrices, and assumes that they are in the appropriate format.
#' It calculates the Fisher Information Matrix (FIM) from the given variational matrix,
#' that is, calculates minus the expectation of the second order derivative matrix of the log-likelihood
#' at the given times with respect to the model parameters.
#' The times are given in the variational matrix, which contains first and second order derivatives of
#' the output with respect to certain parameters, at certain times, for a particular individual.
#'
#' @param df        Data frame of first and second order variations as generated by \code{\link{calcVariations}}
#'   or \code{\link{calcVariationsFim}}, that is, with columns 't', 'i', 'y', 'dy_d<v1>' and
#'   'd2y_d<v1>_<v2>', where variables v1 and v2 are replaced by names.
#'   This should not be normalized. If it is, it is denormalized.
#' @param omega     Symmetric, positive definite, numeric variance-covariance matrix of individual parameters.
#'   Columns and rows should be named, and those names should be the same.
#' @param sigma     Symmetric, positive definite, numeric variance-covariance matrix of residual parameters.
#'   Columns and rows should be named, and those names should be the same.
#' @param vartheta  Vector of names of population parameters for which derivatives are to be calculated.
#'   Should be a subset of the names appearing in the variational matrix, or \code{NULL} for none.
#' @param varomega  Symmetric boolean matrix of the same dimension as \code{omega}, with the same row and
#'   column names.
#'   The value \code{TRUE} means the derivative of the corresponding element should be included in the FIM,
#'   \code{FALSE} means not.
#' @param varsigma  Symmetric boolean matrix of the same dimension as \code{sigma}, with the same row and
#'   column names.
#'   The value \code{TRUE} means the derivative of the corresponding element should be included in the FIM,
#'   \code{FALSE} means not.
#' @param caller    String containing the signature of the calling function, for use in error messages.
#'   A typical value would be 'FisherInfo::calcFimFromMatrix'.
#'
#' @return  The FIM, as a square matrix with element (i,j) equal to -E(d^2 L/d<vi><vj>), where L is the
#'   log likelihood and <vi> is the i-th parameter.
#'   The row and column names of the matrix are those of the variable thetas, omegas and sigmas,
#'   as specified by \code{vartheta}, \code{varomega} and \code{varsigma}.
#'   The FIM is not normalized and the 'thetaNormalized', 'omegaNormalized' and 'sigmaNormalized' attributes
#'   of this matrix are set to named vectors of \code{FALSE}, where the names are the names of the thetas,
#'   and the column names of omegas and sigmas as specified by \code{df}, \code{omega} and \code{sigma},
#'   respectively.
#'   The attributes 'theta', 'omega' and 'sigma' store the values of all parameters as named vectors.
#'   They form a superset of the variable parameters.
#'   The 'type' attribute is set to 'FIM'
#'   The function displays an error and returns \code{NULL} if the input is incorrectly formatted or
#'   the variance matrix cannot be inverted.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
doCalcFim <- function(df, omega, sigma, vartheta, varomega, varsigma, caller) {
  # De-normalize df. This also checks its validity:
  df <- normalizeVariations(df = df, parNorm = FALSE, outNorm = FALSE)
  if (is.null(df)) return(NULL)
  # Identify columns in df of different types of derivatives:
  nam <- checkSetDerivCols(df = df, omega = omega, sigma = sigma, vartheta = vartheta, caller = caller)
  if (is.null(nam)) return(NULL)
  
  # Determine for which omega and sigma elements we need to compute derivatives of the loglikelihood:
  omCases <- createOmegaSigmaCases(varomega, TRUE)          # Data frame with rows, columns and names of omega parameters
  siCases <- createOmegaSigmaCases(varsigma, FALSE)         # Data frame with rows, columns and names of sigma parameters
  nthetas <- length(nam[["theta"]])                         # Nr of fixed-effect parameters (thetas)
  fimdim <- nthetas + nrow(omCases) + nrow(siCases)         # Nr of rows (and cols) of FIM
  fimnam <- c(vartheta, omCases[, "nam"], siCases[, "nam"]) # Row (and col) names of FIM, in order theta, omega, sigma
  
  # Calculate variance matrix Vi for subject i:
  df    <- df[order(df[, "t"], df[, "i"]), ]                # Order by time ("t") and then by output element ("i").
  tjs   <- sort(unique(df[, "t"]))                          # Vector of time points t_j 
  nui   <- as.matrix(df[nam[["eta"]]])                      # Matrix of derivatives nu_i = doutput/deta_i, where i indexes the IIV parameter (eta). Dimension nOutput*nTimes x nEta, where nOutput = nr of output variables, nTimes = nr of time points, nEta = nr of IIV parameters (etas).
  Vi    <- nui %*% omega %*% t(nui)                         # Omega part of variance matrix; dimension nOutput*nTimes x nOutput*nTimes
  nOutput <- length(unique(df[, "i"]))
  for (j in seq_along(tjs)) {
    rhoij <- as.matrix(df[df[, "t"] == tjs[[j]], nam[["eps"]]])  # Matrix of derivatives rho_ij = doutput/deps_ij, where i indexes the residual parameter (eps) and j the time point t_j. Dimension nOutput x nEps, where nEps = nr of residual parameters (eps's).
    rge <- ((j-1)*nOutput + 1) : (j*nOutput)                     # Index range of the submatrix of Vi corresponding to the j-th time point.
    Vi[rge, rge] <- Vi[rge, rge] + rhoij %*% sigma %*% t(rhoij)  # Add Sigma part of variance matrix per time point to Vi
  }
  
  # Compute inverse of Vi:
  scaleVi  <- FALSE               # Flag whether Vi is rescaled. Optionally Vi is rescaled by 1/output to facilitate calculation of inverse (this rescaling removes the size of the output y from the matrix Vi):
  scaleMat <- NULL                # Rescaling matrix (if Vi is rescaled).
  Vinv <- try(solve(Vi), silent = TRUE)
  if(inherits(Vinv, "try-error")) {
    err <- FALSE
    scaleVi  <- TRUE
    scaleMat <- diag(1/df[, "y"], nrow = length(df[, "y"]))
    if (any(is.infinite(scaleMat))) {
      err <- TRUE
    } else {
      Vinv <- try(solve(scaleMat %*% Vi %*% scaleMat), silent = TRUE)
    }
    if(err || inherits(Vinv, "try-error") || any(is.na(Vinv))) {
      processErrors(paste0(caller, ": Variance matrix cannot be inverted, not even after rescaling.\nExiting.\n"))
      return(NULL)
    }
    processWarns(paste0(caller, ": Variance matrix had to be rescaled to be inverted. No need to worry, I can continue.\n"))
  }
  
  # Compute inverse(Vi) * dVi/dpar for par in theta, omega or sigma:
  VinvVth  <- lapply(vartheta, function(th) {
    thEta <- grep(paste0("_d", th, "_|_", th, "$"), nam[["thetaEta"]], value = TRUE) # Names of applicable theta-eta 2nd derivs
    thEps <- grep(paste0("_d", th, "_|_", th, "$"), nam[["thetaEps"]], value = TRUE) # Names of applicable theta-eps 2nd derivs
    dnuidth <- as.matrix(df[thEta])  # Matrix of derivatives dnu_i/dtheta = d^2output/(dtheta deta_i), for the theta given by th; nu_i is as above. Dimension nOutput*nTimes x nEta
    Vth <- dnuidth %*% omega %*% t(nui) + nui %*% omega %*% t(dnuidth) # Omega part of dVi/dtheta, for the theta given by th; dimension nOutput*nTimes x nOutput*nTimes
    for (j in seq_along(tjs)) { # j <- 1
      rhoij         <- as.matrix(df[df[, "t"] == tjs[[j]], nam[["eps"]]])  # Matrix of derivatives rho_ij = doutput/deps_ij, as above. Dimension nOutput x nEps
      drhoijdth     <- as.matrix(df[df[, "t"] == tjs[[j]], thEps])         # Matrix of derivatives drho_ij/dtheta = d^2output/(dtheta deps_ij), for the theta given by th. Dimension nOutput x nEps
      rge           <- ((j-1)*nOutput + 1) : (j*nOutput)                   # Index range of the submatrix of Vi corresponding to the j-th time point, as above.
      Vth[rge, rge] <- Vth[rge, rge] + drhoijdth %*% sigma %*% t(rhoij) + rhoij %*% sigma %*% t(drhoijdth)
    } # add Sigma part of derivate dV/dth per time point to Vth
    return(Vinv %*% if(scaleVi) scaleMat %*% Vth %*% scaleMat else Vth)
  })  # Contains a list with elements inv(V) * dV/dtheta for all thetas in vartheta. May be empty.
  VinvVom <- plyr::alply(omCases, .margins = 1, .fun = function(case) {
    rw <- case[1, "row"]           # Row index of the Omega element described by "case"
    cl <- case[1, "col"]           # Column index of the Omega element described by "case"
    xi = if(rw == cl) 0.5 else 1   # Factor 1/2 if derivative is to be taken wrt diagonal elt of Omega
    Vom <- nui[, rw] %*% t(nui[, cl]) + nui[, cl] %*% t(nui[, rw])
    return(xi * Vinv %*% if(scaleVi) scaleMat %*% Vom %*% scaleMat else Vom)
  })  # Contains a list with elements inv(V) * dV/domega for all omegas in varomega. May be empty.
  VinvVsi <- plyr::alply(siCases, .margins = 1, .fun = function(case) {
    rw <- case[1, "row"]           # Row index of the Sigma element described by "case"
    cl <- case[1, "col"]           # Column index of the Sigma element described by "case"
    xi = if(rw == cl) 0.5 else 1   # Factor 1/2 if derivative is to be taken wrt diagonal elt of Sigma
    Vsi <- 0 * Vi  # zero matrix with dimension of Vi
    for (j in seq_along(tjs)) {
      rhoij <- as.matrix(df[df[, "t"] == tjs[[j]], nam[["eps"]]])
      rge <- ((j-1)*nOutput + 1) : (j*nOutput)
      Vsi[rge, rge] <- xi * (rhoij[, rw] %*% t(rhoij[, cl]) + rhoij[, cl] %*% t(rhoij[, rw]))
    }
    return(Vinv %*% if(scaleVi) scaleMat %*% Vsi %*% scaleMat else Vsi)
  })  # Contains a list with elements inv(V) * dV/dsigma for all sigmas in varsigma. May be empty.
  VinvVpar <- c(VinvVth, VinvVom, VinvVsi)  # List of matrices inv(V) * dV/dpar where par is all variable theta, omega, sigma.

  # Compute t(doutput/dtheta) Vinv doutput/dtheta:  
  mu <- as.matrix(df[nam[["theta"]]])   # Matrix mu = doutput/dtheta. Dimension nOutput*nTimes x nthetas
  if(scaleVi) mu <- scaleMat %*% mu
  muVinvmu <- t(mu) %*% Vinv %*% mu     # The nthetas x nthetas matrix with elements (doutput/dtheta)^t * inv(V) * doutput/dtheta for all thetas in vartheta. May be empty.
  
  # Create empty FIM:
  fim <- matrix(0, nrow = fimdim, ncol = fimdim, dimnames = list(fimnam, fimnam))
  # Fill FIM:
  for (i in 1:fimdim) {  # i <- 1
    for (j in 1:i) {     # j <- 1
      fim[i, j] <- 0.5 * traceprod(t(VinvVpar[[i]]), VinvVpar[[j]]) + if (i <= nthetas) muVinvmu[i, j] else 0 
      fim[j, i] <- fim[i, j]
    }
  }
  attr(fim, typeAttr)      <- varTypeFim
  attr(fim, symbAttr)      <- attr(df, symbAttr)
  theta                    <- attr(df, thetaAttr)
  attr(fim, thetaAttr)     <- theta
  attr(fim, omegaAttr)     <- omega
  attr(fim, sigmaAttr)     <- sigma
  thetaNorm <- rep(FALSE, length(theta))
  names(thetaNorm) <- names(theta)
  omegaNorm <- rep(FALSE, nrow(omega))
  names(omegaNorm) <- colnames(omega)
  sigmaNorm <- rep(FALSE, nrow(sigma))
  names(sigmaNorm) <- colnames(sigma)
  attr(fim, thetaNormAttr) <- thetaNorm
  attr(fim, omegaNormAttr) <- omegaNorm
  attr(fim, sigmaNormAttr) <- sigmaNorm
  class(fim)               <- c('keepattr', class(fim))   # This keeps attributes when subsetting
  return(fim)
}


#------------------------------------------ isValidFim ------------------------------------------
#' Check validity of a FIM
#'
#' Checks whether the provided input is a valid FIM,
#' that is, whether it is in the format as produced by \code{\link{calcFimFromMatrix}} or \code{\link{calcFimFromModel}}.
#' In particular, it is checked whether the input is a matrix of numeric values, with attributes 'type', 'symbolic',
#' 'omega', 'sigma', 'thetaNormalized', 'omegaNormalized' and 'sigmaNormalized', where 'type' is set to 'FIM' and
#' 'symbolic' is a boolean, and whether the attribute names are consistent.
#' It is not checked whether attribute 'theta' is present, because this parameter may be NULL in which case it is
#' not listed.
#'
#' @param df     Any object.
#' @param nmdf   Name (string) by which this object should be referred in any error messages. By default "input 'df'".
#'
#' @return \code{NULL} if \code{df} is a valid FIM,
#'   and an error message or list of error messages if not.
#'
#'
#' @author Martijn van Noort
#' 
#' @noRd 
isValidFim <- function(df, nmdf = "input 'df'") {
  out <- NULL
  if(!"keepattr" %in% class(df)) {
    return(paste0(nmdf, " is not of the 'keepattr' class."))
  }
  if(!all(c(typeAttr, symbAttr, omegaAttr, sigmaAttr, thetaNormAttr, omegaNormAttr, sigmaNormAttr)
          %in% names(attributes(df))) || attr(df, typeAttr) != varTypeFim || !is.logical(attr(df, symbAttr))) {
    return(paste0(nmdf, " does not have the correct attributes of a FIM."))
  }
  if (!is.matrix(df)) {
    return(paste0(nmdf, " should be a matrix."))
  }
  if (!is.numeric(df)) {
    out <- c(out, list(paste0(nmdf, " contains non-numeric values.")))
  }
  if(nrow(df) != ncol(df)) {
    out <- c(out, list(paste0(nmdf, " is not square.")))
  }
  if(nrow(df) != length(colnames(df)) || nrow(df) != length(row.names(df))) {
    out <- c(out, list(paste0(nmdf, " has unnamed rows or columns.")))
  }
  if(nrow(df) > 0 && any(colnames(df) != row.names(df))) {
    out <- c(out, list(paste0(nmdf, " has non-matching row and column names.")))
  }
  att1 <- attr(df, thetaNormAttr)
  att2 <- attr(df, thetaAttr)
  if((length(att1) > 0 || length(att2) > 0) && any(names(att1) != names(att2))) {
    out <- c(out, list(paste0(nmdf, " attributes 'thetaAttr' and 'thetaNormAttr' are inconsistently named.")))
  }
  att1 <- attr(df, omegaNormAttr)
  att2 <- attr(df, omegaAttr)
  if((length(att1) > 0 || ncol(att2) > 0) && any(names(att1) != colnames(att2))) {
    out <- c(out, list(paste0(nmdf, " attributes 'omegaAttr' and 'omegaNormAttr' are inconsistently named.")))
  }
  att1 <- attr(df, sigmaNormAttr)
  att2 <- attr(df, sigmaAttr)
  if((length(att1) > 0 || ncol(att2) > 0) && any(names(att1) != colnames(att2))) {
    out <- c(out, list(paste0(nmdf, " attributes 'sigmaAttr' and 'sigmaNormAttr' are inconsistently named.")))
  }
  # Check names of FIM:
  nams <- unlist(lapply(c(omegaAttr, sigmaAttr), function(matAttr) {
    varmat <- apply(attr(df, matAttr), 1:2, function(x) TRUE)  # Create matrix like omega or sigma, with all elements set to TRUE
    createOmegaSigmaCases(varmat = varmat, isOmega = TRUE)[, "nam"]
  }))
  nams <- c(names(attr(df, thetaAttr)), nams)  # These are all potential names
  if(!all(row.names(df) %in% nams)) {
    out <- c(out, list(paste0(nmdf, " has invalid row and column names.")))
  }
  return(out)
}


#--------------------------------------------------------------------------------------------------------
#------------------------------------------ Exported functions ------------------------------------------
#--------------------------------------------------------------------------------------------------------


#------------------------------------------ calcFimFromMatrix ------------------------------------------
#' Calculate Fisher Information Matrix
#'
#' Calculates the Fisher Information Matrix (FIM) from the given variational matrix,
#' that is, calculates minus the expectation of the second order derivative matrix of the log-likelihood
#' at the given times with respect to the model parameters.
#' The times are given in the variational matrix, which contains first and second order derivatives of
#' the output with respect to certain parameters, at certain times, for a particular individual.
#'
#' @param df        Data frame of first and second order variations as generated by \code{\link{calcVariations}}
#'   or \code{\link{calcVariationsFim}}, that is, with columns 't', 'i', 'y', 'dy_d<v1>' and
#'   'd2y_d<v1>_<v2>', where variables v1 and v2 are replaced by names.
#' @param omega     Variance-covariance matrix of individual parameters, as numeric matrix.
#'   Should be symmetric and positive definite.
#'   Columns and/or rows should be named. If both, the names should be the same.
#' @param sigma     Variance-covariance matrix of residual parameters, as numeric matrix.
#'   Should be symmetric and positive definite.
#'   Columns and/or rows should be named. If both, the names should be the same.
#' @param vartheta  Vector of names of population parameters for which derivatives are to be calculated.
#'   Should be a subset of the names appearing in the variational matrix, or \code{NULL} for none.
#' @param varomega  Matrix (not necessarily named) of booleans of the same dimension as \code{omega}.
#'   The value \code{TRUE} means the derivative of the corresponding element should be included in the FIM,
#'   \code{FALSE} means not.
#'   The matrix should be symmetric.
#'   Instead of a matrix, may also provide a single boolean or \code{NULL}, where \code{TRUE} stands for a
#'   matrix of all \code{TRUE} (also the off-diagonal elements), \code{FALSE} for a matrix of all \code{FALSE},
#'   and \code{NULL} (default) is to set all nonzero elements of \code{omega} to TRUE, and the others to \code{FALSE}.
#' @param varsigma  Matrix (not necessarily named) of booleans of the same dimension as \code{sigma}.
#'   The value \code{TRUE} means the derivative of the corresponding element should be included in the FIM,
#'   \code{FALSE} means not.
#'   The matrix should be symmetric.
#'   Instead of a matrix, may also provide a single boolean or \code{NULL}, where \code{TRUE} stands for a
#'   matrix of all \code{TRUE} (also the off-diagonal elements), \code{FALSE} for a matrix of all \code{FALSE},
#'   and \code{NULL} (default) is to set all nonzero elements of \code{sigma} to TRUE, and the others to \code{FALSE}.
#'
#' @return  The FIM, as a square matrix with element (i,j) equal to -E(d^2 L/d<vi><vj>), where L is the
#'   log likelihood and <vi> is the i-th parameter.
#'   The row and column names of the matrix are those of the variable thetas, omegas and sigmas,
#'   as specified by \code{vartheta}, \code{varomega} and \code{varsigma}.
#'   The FIM is not normalized and the 'thetaNormalized', 'omegaNormalized' and 'sigmaNormalized' attributes
#'   of this matrix are set to named vectors of \code{FALSE}, where the names are the names of the thetas,
#'   and the column names of omegas and sigmas as specified by \code{df}, \code{omega} and \code{sigma},
#'   respectively.
#'   The attributes 'theta', 'omega' and 'sigma' store the values of all parameters as named vectors.
#'   They form a superset of the variable parameters.
#'   The 'type' attribute is set to 'FIM'
#'   The function displays an error and returns \code{NULL} if the input is incorrectly formatted or
#'   the variance matrix cannot be inverted.
#'
#' @export
#' 
#' @family calculation functions
#'
#' @author Martijn van Noort
calcFimFromMatrix <- function(df, omega, sigma, vartheta, varomega = NULL, varsigma = NULL) {
  # Checks and formatting:
  out <- checkSetCovMat(omega, sigma, varomega, varsigma, caller = 'FisherInfo::calcFimFromMatrix')
  if (is.null(out)) return(NULL)
  omega <- out[[1]]
  sigma <- out[[2]]
  varomega <- out[[3]]
  varsigma <- out[[4]]

  doCalcFim(df = df, omega = omega, sigma = sigma, vartheta = vartheta, varomega = varomega, varsigma = varsigma,
            caller = 'FisherInfo::calcFimFromMatrix')
}


#------------------------------------------ calcFimFromModel ------------------------------------------
#' Calculate Fisher Information Matrix
#'
#' Calculates the Fisher Information Matrix (FIM) from the given model,
#' that is, calculates minus the expectation of the second order derivative matrix of the log-likelihood
#' at the given times with respect to the model parameters.
#' These parameters are listed in \code{vartheta}, \code{varomega} and \code{varsigma}.
#'
#' @param model     Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
#'   This function should return the numeric vector dy/dt.
#' @param p         Function(theta, eta) of population and individual parameters specifying the model parameters.
#'   This function should return a named numeric vector.
#' @param init      Function(p) of parameters specifying the initial state of the model.
#' @param output    Function(y, p, eps) of state variables, parameters and residual errors specifying the model outputs.
#'   This function should return a numeric vector.
#' @param times     Numeric vector of times where variations are to be evaluated.
#' @param theta     Population parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param omega     Variance-covariance matrix of individual parameters, as numeric matrix.
#'   Should be symmetric and positive definite.
#'   Columns and/or rows should be named. If both, the names should be the same.
#' @param sigma     Variance-covariance matrix of residual parameters, as numeric matrix.
#'   Should be symmetric and positive definite.
#'   Columns and/or rows should be named. If both, the names should be the same.
#' @param vartheta  Vector of names of population parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(theta)}, or \code{NULL} for none.
#'   By default equal to \code{names(theta)}.
#' @param varomega  Matrix (not necessarily named) of booleans of the same dimension as \code{omega}.
#'   The value \code{TRUE} means the derivative of the corresponding element should be included in the FIM,
#'   \code{FALSE} means not.
#'   The matrix should be symmetric.
#'   Instead of a matrix, may also provide a single boolean or \code{NULL}, where \code{TRUE} stands for a
#'   matrix of all \code{TRUE} (also the off-diagonal elements), \code{FALSE} for a matrix of all \code{FALSE},
#'   and \code{NULL} (default) is to set all nonzero elements of \code{omega} to TRUE, and the others to \code{FALSE}.
#' @param varsigma  Matrix (not necessarily named) of booleans of the same dimension as \code{sigma}.
#'   The value \code{TRUE} means the derivative of the corresponding element should be included in the FIM,
#'   \code{FALSE} means not.
#'   The matrix should be symmetric.
#'   Instead of a matrix, may also provide a single boolean or \code{NULL}, where \code{TRUE} stands for a
#'   matrix of all \code{TRUE} (also the off-diagonal elements), \code{FALSE} for a matrix of all \code{FALSE},
#'   and \code{NULL} (default) is to set all nonzero elements of \code{sigma} to TRUE, and the others to \code{FALSE}.
#' @param symbolic \code{TRUE} (default) if derivatives are to be computed symbolically, \code{FALSE} if numerically.
#'   See \code{\link{calcVariations}} for notes on using symbolic derivation.
#' @param chkModel \code{TRUE} (default) if it has to be checked whether model components (\code{model}, \code{p}, \code{init}, \code{output}) 
#'   are formatted correctly for symbolic manipulation, \code{FALSE} if not.
#' @param ...    Named arguments to be passed to \code{\link[deSolve]{lsoda}}.
#'   Can be used for example to pass events or error tolerances.
#'
#' @return  The FIM, as a square matrix with element (i,j) equal to -E(d^2 L/d<vi><vj>), where L is the
#'   log likelihood and <vi> is the i-th parameter.
#'   The row and column names of the matrix are those of the variable thetas, omegas and sigmas,
#'   as specified by \code{vartheta}, \code{varomega} and \code{varsigma}.
#'   The FIM is not normalized and the 'thetaNormalized', 'omegaNormalized' and 'sigmaNormalized' attributes
#'   of this matrix are set to named vectors of \code{FALSE}, where the names are the names of the thetas,
#'   and the column names of omegas and sigmas as specified by \code{df}, \code{omega} and \code{sigma},
#'   respectively.
#'   The attributes 'theta', 'omega' and 'sigma' store the values of all parameters as named vectors.
#'   They form a superset of the variable parameters.
#'   The 'type' attribute is set to 'FIM', and the 'symbolic' one to the value of \code{symbolic}.
#'   The function displays an error and returns \code{NULL} if the input is incorrectly formatted or
#'   the variance matrix cannot be inverted.
#'
#' @note It is recommended to use symbolic derivation (\code{symbolic==TRUE}) rather than numerical derivation,
#'   as it is more accurate, especially for derivatives that are formally zero.
#'   See \code{\link{calcVariations}} for details.
#'   
#' @export
#'
#' @family calculation functions
#' 
#' @author Martijn van Noort
calcFimFromModel <- function(model, p, init, output, times, theta, omega, sigma, vartheta = names(theta),
                             varomega = NULL, varsigma = NULL, symbolic = TRUE, chkModel = TRUE, ...) {
  # Checks and formatting:
  out <- checkSetCovMat(omega, sigma, varomega, varsigma, caller = 'FisherInfo::calcFimFromModel')
  if (is.null(out)) return(NULL)
  omega <- out[[1]]
  sigma <- out[[2]]
  varomega <- out[[3]]
  varsigma <- out[[4]]
  nmeta <- colnames(omega)
  nmeps <- colnames(sigma)
  df <- calcVariationsFim(model = model, p = p, init = init, output = output, times = times, theta = theta,
                          nmeta = nmeta, nmeps = nmeps, vartheta = vartheta, vareta = nmeta, vareps = nmeps,
                          symbolic = symbolic, chkModel = chkModel, ...)
  doCalcFim(df = df, omega = omega, sigma = sigma, vartheta = vartheta, varomega = varomega, varsigma = varsigma,
            caller = 'FisherInfo::calcFimFromModel')
}


#------------------------------------------ getAllParsFim ------------------------------------------
#' Get all parameter values of a FIM
#'
#' Gets the values of the parameters used to create a given FIM.
#'
#' @param fim    The FIM.
#'
#' @return A named list of three named numeric vectors containing the theta's, omega's and sigma's used to
#'   create the FIM.
#'   The list elements are named "theta", "omega" and "sigma", respectively.
#'   The vector elements have the names of the parameters, where off-diagonal elements of omega and sigma are
#'   named "(<par1>,<par2>)".
#'   The function displays an error and returns \code{NULL} if \code{fim} is not a valid FIM.
#'
#' @export
#'
#' @family retrieval functions
#' 
#' @author Martijn van Noort
getAllParsFim <- function(fim) {
  # Check input validity:
  valid <- isValidFim(fim, nmdf = "input 'fim'")
  if (!is.null(valid)) {
    processErrors(paste0("FisherInfo::getAllParsFim:\n", paste0(paste0("\t", unlist(valid), "\n"), collapse = ""), "Exiting.\n"))
    return(NULL)
  }
  # Get omega and sigma values 
  vals <- lapply(c(omegaAttr, sigmaAttr), function(matAttr) {
    mat <- attr(fim, matAttr)
    varmat <- apply(mat, 1:2, function(x) TRUE)  # Create matrix like omega or sigma, with all elements set to TRUE
    allValues <- createOmegaSigmaCases(varmat = varmat, isOmega = TRUE)
    allValues[, "value"] <- apply(allValues, 1, function(row) mat[as.numeric(row["row"]), as.numeric(row["col"])])
    out <- allValues[, "value"]
    names(out) <- allValues[, "nam"]
    return(out)
  })
  vals <- c(list(attr(fim, thetaAttr)), vals)
  names(vals) <- c("theta", "omega", "sigma")
  return(vals)
}


#------------------------------------------ getParsFim ------------------------------------------
#' Get the column parameter values of a FIM
#'
#' Gets the values of the parameters corresponding to the columns of the FIM.
#' This is a subset of \code{getAllParsFim} corresponding to the columns present in the FIM.
#'
#' @param fim    The FIM.
#'
#' @return A named list of three named numeric vectors containing the theta's, omega's and sigma's used to
#'   create the FIM, and present as columns of the FIM, in the order of appearance in the FIM.
#'   The list elements are named "theta", "omega" and "sigma", respectively.
#'   The vector elements have the names of the parameters, where off-diagonal elements of omega and sigma are
#'   named "(<par1>,<par2>)".
#'   The function displays an error and returns \code{NULL} if \code{fim} is not a valid FIM.
#'
#' @export
#' 
#' @family retrieval functions
#'
#' @author Martijn van Noort
getParsFim <- function(fim) {
  out <- getAllParsFim(fim)
  if(is.null(out)) return(NULL)
  lapply(out, function(vec) vec[colnames(fim)[colnames(fim) %in% names(vec)]])
}


#------------------------------------------ getParVecFim ------------------------------------------
#' Get the column parameter values of a FIM as a single vector
#'
#' Gets the values of the parameters corresponding to the columns of the FIM.
#' This is a subset of \code{getAllParsFim} corresponding to the columns present in the FIM.
#'
#' @param fim    The FIM.
#'
#' @return A named numeric vector containing the theta's, omega's and sigma's used to
#'   create the FIM, and present as columns of the FIM, in the order of appearance in the FIM.
#'   The names are the names of the parameters, where off-diagonal elements of omega and sigma are
#'   named "(<par1>,<par2>)".
#'   The function displays an error and returns \code{NULL} if \code{fim} is not a valid FIM.
#'
#' @export
#' 
#' @family retrieval functions
#'
#' @author Martijn van Noort
getParVecFim <- function(fim) {
  out <- getAllParsFim(fim)
  if(is.null(out)) return(NULL)
  names(out) <- NULL
  out <- unlist(out, recursive = FALSE, use.names = TRUE)
  if(!all(colnames(fim) %in% names(out))) {
    # Unknown columns in FIM. This should not happen:
    processErrors(paste0("FisherInfo::getParVecFim: no parameter values for columns ",
                         paste0(setdiff(colnames(fim), names(out)), ", "), ".\nExiting.\n"))
    return(NULL)
  }
  out[colnames(fim)]
}


#------------------------------------------ normalizeFim ------------------------------------------
#' Normalize Fisher Information Matrix
#'
#' Normalizes the given Fisher Information Matrix (FIM) with respect to the given parameters.
#' Effectively this means that the derivatives with respect to these parameters are taken relative to their size, i.e.,
#' they are taken with respect to the logarithms of these parameters.
#' This should be used for those parameters for which relative changes are relevant, rather than absolute changes.
#' Normalization is accomplished by left and right multiplication of the given FIM by a diagonal matrix of normalizing factors.
#' The normalizing factors are the values of the corresponding parameters, and are also given as input.
#' If the given FIM was already normalized (as indicated by the "thetaNormalized", "omegaNormalized" and
#' "sigmaNormalized" attributes) then this normalization is undone.
#'
#' @param fim        A matrix with row and column names containing the FIM,
#'   as produced by \code{\link{calcFimFromMatrix}} or \code{\link{calcFimFromModel}}.
#' @param thetaNorm  Vector (not necessarily named) of booleans of the same dimension as the theta vector
#'   used to create the FIM (so not just the variable theta's).
#'   The value \code{TRUE} means the corresponding structural parameter should be normalized in the FIM,
#'   \code{FALSE} means not.
#'   Instead of a vector, may also provide a single boolean, where \code{TRUE} (default) stands for a vector of
#'   all \code{TRUE}, and \code{FALSE} for a vector of all \code{FALSE}.
#'   May also provide a character vector with the names of the parameters to be normalized.
#'   Invalid names are ignored, with a warning.
#' @param omegaNorm  Vector (not necessarily named) of booleans of the same dimension as the omega matrix
#'   used to create the FIM (so not just the variable omega's).
#'   The value \code{TRUE} means the corresponding individual parameter should be normalized in the FIM,
#'   \code{FALSE} means not.
#'   Instead of a vector, may also provide a single boolean, where \code{TRUE} (default) stands for a vector of
#'   all \code{TRUE}, and \code{FALSE} for a vector of all \code{FALSE}.
#'   May also provide a character vector with the names of the parameters to be normalized.
#'   Invalid names are ignored, with a warning.
#' @param sigmaNorm  Vector (not necessarily named) of booleans of the same dimension as the sigma matrix
#'   used to create the FIM (so not just the variable sigma's).
#'   The value \code{TRUE} means the corresponding residual parameter should be normalized in the FIM,
#'   \code{FALSE} means not.
#'   Instead of a vector, may also provide a single boolean, where \code{TRUE} (default) stands for a vector of
#'   all \code{TRUE}, and \code{FALSE} for a vector of all \code{FALSE}.
#'   May also provide a character vector with the names of the parameters to be normalized.
#'   Invalid names are ignored, with a warning.
#' @param checkPos   If \code{TRUE} (default), it is checked that parameters to be normalized in the FIM
#'   are strictly positive. If \code{FALSE}, this is not checked.
#'
#' @return   A matrix in same format as \code{fim}, where the values are normalized as specified, and the
#'   attributes are adapted accordingly.
#'   The function displays an error and returns \code{NULL} if the input is not in the specified format.
#'
#' @details This can be used for parameters for which relative changes in value are the relevant quantity,
#'   rather than absolute changes, e.g., clearances, rates, volumes, and additive errors.
#'   Absolute changes in such parameters are not relevant because they depend on the chosen units.
#'   Relative changes do not suffer from this problem.
#'   It should typically also be used for parameters that already have a relative size,
#'   such as IIV variances (in case of multiplicative IIV effects, i.e., as exp(eta) factors),
#'   proportional errors, powers ("gamma"), Imax, etc.
#'   The reason is that these parameters may also be scaled by arbitrary multiplication factors
#'   (e.g., a proportional error can be modelled as a fraction or as a percentage).
#'
#' @export
#' 
#' @family result modifiers
#' 
#' @author Martijn van Noort
normalizeFim <- function(fim, thetaNorm = TRUE, omegaNorm = TRUE, sigmaNorm = TRUE, checkPos = TRUE) {
  valid <- isValidFim(fim, nmdf = "input 'fim'")
  if (!is.null(valid)) {
    processErrors(paste0("FisherInfo::normalizeFim:\n", paste0(paste0("\t", unlist(valid), "\n"), collapse = ""), "Exiting.\n"))
    return(NULL)
  }
  theta <- attr(fim, thetaAttr)
  omega <- attr(fim, omegaAttr)
  sigma <- attr(fim, sigmaAttr)
  # Checks and formatting:
  out <- checkSetCovMat(omega, sigma, FALSE, FALSE, caller = 'FisherInfo::normalizeFim')
  if (is.null(out)) return(NULL)
  omega <- out[[1]]
  sigma <- out[[2]]
  # Set thetaNorm to named vector of boolean:
  if (is.logical(thetaNorm)) {
    if (length(thetaNorm) == 1) thetaNorm <- rep(thetaNorm, length(theta))
  } else {
    # thetaNorm treated as character vector:
    if (!all(thetaNorm %in% names(theta))) processWarns("FisherInfo::normalizeFim: 'thetaNorm' contains unknown parameters. They will be ignored.\nContinuing.\n")
    thetaNorm <- names(theta) %in% thetaNorm
  }
  if (length(thetaNorm) != length(theta)) {
    processErrors(paste0("FisherInfo::normalizeFim: 'thetaNorm' has incorrect format.\nExiting.\n"))
    return(NULL)
  }
  names(thetaNorm) <- names(theta)
  # Set omegaNorm to named vector of boolean:
  if (is.logical(omegaNorm)) {
    if (length(omegaNorm) == 1) omegaNorm <- rep(omegaNorm, nrow(omega))
  } else {
    # omegaNorm treated as character vector:
    if (!all(omegaNorm %in% colnames(omega))) processWarns("FisherInfo::normalizeFim: 'omegaNorm' contains unknown parameters. They will be ignored.\nContinuing.\n")
    omegaNorm <- colnames(omega) %in% omegaNorm
  }
  if (length(omegaNorm) != nrow(omega)) {
    processErrors(paste0("FisherInfo::normalizeFim: 'omegaNorm' has incorrect format.\nExiting.\n"))
    return(NULL)
  }
  names(omegaNorm) <- colnames(omega)
  # Set sigmaNorm to named vector of boolean:
  if (is.logical(sigmaNorm)) {
    if (length(sigmaNorm) == 1) sigmaNorm <- rep(sigmaNorm, nrow(sigma))
  } else {
    # sigmaNorm treated as character vector:
    if (!all(sigmaNorm %in% colnames(sigma))) processWarns("FisherInfo::normalizeFim: 'sigmaNorm' contains unknown parameters. They will be ignored.\nContinuing.\n")
    sigmaNorm <- colnames(sigma) %in% sigmaNorm
  }
  if (length(sigmaNorm) != nrow(sigma)) {
    processErrors(paste0("FisherInfo::normalizeFim: 'sigmaNorm' has incorrect format.\nExiting.\n"))
    return(NULL)
  } 
  names(sigmaNorm) <- colnames(sigma)
  # Determine for which omegas and sigmas we need to normalize, and determine the parameter values:
  existThetaNorm <- attr(fim, thetaNormAttr)
  existOmegaNorm <- attr(fim, omegaNormAttr)
  existSigmaNorm <- attr(fim, sigmaNormAttr)
  omCases <- computeNormalizers(omega, omegaNorm, existOmegaNorm)
  siCases <- computeNormalizers(sigma, sigmaNorm, existSigmaNorm)
  # Determine normalizing values for theta:
  thvals <- theta
  thvals[!thetaNorm] <- 1
  # Check all thetas have normalizers (this was checked in computeNormalizers for omega and sigma):
  if(!all(names(existThetaNorm) %in% names(thvals))) {
    processErrors(paste0("FisherInfo::normalizeFim: parameters ",
                         paste0(setdiff(names(thvals), names(existThetaNorm)), collapse = ", "),
                         " not specified in 'theta'. Exiting.\n"))
    return(NULL)
  }
  if(length(existThetaNorm) > 0) {
    # Undo prior normalization:
    fac <- thvals[existThetaNorm]
    fac[fac <= 0] <- 1  # if fac == 0 then prior normalization would have normalized to 0, in which case the value of 'fac' does not matter.
    thvals[existThetaNorm] <- thvals[existThetaNorm]/fac
  }
  # Combine normalizers:
  normVals <- c(thvals, omCases[, "value"], siCases[, "value"])
  names(normVals) <- c(names(thvals), omCases[, "nam"], siCases[, "nam"])
  # Check that all elements of FIM have a normalizing values. This should be the case as it was already checked above:
  if (length(setdiff(colnames(fim), names(normVals))) > 0) {
    processErrors(paste0("FisherInfo::normalizeFim: no normalizing value for ",
               paste0(setdiff(colnames(fim), names(normVals)), collapse = ", "), ". Exiting.\n"))
    return(NULL)
  }
  normVals <- normVals[colnames(fim)]   # Restrict to elements occurring in FIM, in the same order as in the FIM
  # Optionally check that normalizing values are positive:
  if (checkPos & any(normVals <= 0)) {
    processErrors("FisherInfo::normalizeFim: normalizing value is negative. Exiting.\n")
    return(NULL)
  }
  # Perform normalization:
  normalizer <- diag(normVals, nrow = nrow(fim))
  row.names(normalizer) <- colnames(normalizer) <- colnames(fim)
  nfim <- normalizer %*% fim %*% normalizer  # This loses attributes, so they have to be set again:
  attributes(nfim) <- attributes(fim)
  class(nfim) <- class(fim)  # Not set by "attributes"
  attr(nfim, thetaNormAttr) <- thetaNorm
  attr(nfim, omegaNormAttr) <- omegaNorm
  attr(nfim, sigmaNormAttr) <- sigmaNorm
  return(nfim)
}


#------------------------------------------ isNormalizedFim ------------------------------------------
#' Checks whether a FIM is normalized
#'
#' Returns whether the given FIM is normalized, as recorded in its normalization indicators.
#'
#' @param fim        A matrix with row and column names containing the FIM,
#'   as produced by \code{\link{calcFimFromMatrix}} or \code{\link{calcFimFromModel}}.
#'
#' @return   a named boolean vector, stating for each parameter of the FIM whether it is
#'   normalized (\code{TRUE}) or not (\code{FALSE}).
#'   All theta, omega and sigma parameters are listed, not just the ones present in the FIM.
#'   Returns \code{NULL} and prints an error if there was no normalization indicator.
#'
#' @export
#' 
#' @family checkers
#' 
#' @author Martijn van Noort
isNormalizedFim <- function(fim) {
  valid <- isValidFim(fim, nmdf = "input 'fim'")
  if (!is.null(valid)) {
    processErrors(paste0("FisherInfo::isNormalizedFim:\n", paste0(paste0("\t", unlist(valid), "\n"), collapse = ""), "Exiting.\n"))
    return(NULL)
  }
  out <- c(attr(fim, thetaNormAttr), attr(fim, omegaNormAttr), attr(fim, sigmaNormAttr))
  if (is.null(out)) {
    processErrors("FisherInfo::isNormalizedFim: normalization not recorded.\nExiting.\n")
    return(FALSE)
  }
  return(out)
}


#------------------------------------------ subFim ------------------------------------------
#' Subset a FIM matrix
#'
#' \code{subFim(fim, inds)} works just like \code{fim[inds, inds]}.
#' Use \code{subFim} instead of normal matrix subsetting functions to preserve attributes.
#' 
#' @param fim  A matrix with row and column names containing the FIM,
#'   as produced by \code{\link{calcFimFromMatrix}} or \code{\link{calcFimFromModel}}.
#' @param inds Vector of indices of elements to keep. They may be provided as numerals or names.
#'
#' @return   The submatrix \code{fim[inds, inds]}.
#'   Returns \code{NULL} and prints an error if the input \code{fim} is not a valid FIM.
#'   Returns the input \code{fim} and prints an error if the result of subsetting is not a valid FIM.
#'   This should not be possible.
#'
#' @export
#' 
#' @family retrieval functions
#'
#' @author Martijn van Noort
subFim <- function(fim, inds) {
  valid <- isValidFim(fim, nmdf = "input 'fim'")
  if (!is.null(valid)) {
    processErrors(paste0("FisherInfo::subFim:\n", paste0(paste0("\t", unlist(valid), "\n"), collapse = ""), "Exiting.\n"))
    return(NULL)
  }
  fim[inds, inds]
}


#------------------------------------------ fimIdent ------------------------------------------
#' Determine identifiability using the Fisher Information Matrix
#'
#' Determines identifiability for a given Fisher Information Matrix (FIM).
#' To this end, principal directions are determined in the parameter space, and ordered by curvature.
#' The least curved direction is the direction in which the Objective Function Value (OFV) increases
#' the least with a unit change in parameters. 
#' A direction is considered non-identifiable if its curvature is not larger than a given threshold.
#' 
#' Optionally, the FIM is normalized with respect to certain parameters.
#' Effectively this means that derivatives with respect to these parameters are taken relative to their size, i.e.,
#' they are taken with respect to the logarithms of these parameters.
#' This should be used for those parameters for which relative changes are relevant, rather than absolute changes.
#' For example, one could normalize those parameters that have a unit, but not the ones that have a natural scale
#' (such as a parameter that models a fraction or an exponent).
#' Normalization is accomplished by left and right multiplication of the given FIM by a diagonal matrix of
#' normalizing factors.
#' The normalizing factors are the values of the corresponding parameters, and are also given as input.
#'
#' @param fim             A matrix with row and column names containing the FIM,
#'   as produced by \code{\link{calcFimFromMatrix}} or \code{\link{calcFimFromModel}}.
#' @param curvature       Threshold curvature. Parameter directions with curvatures lower than this value are
#'   considered non-identifiable.
#' @param thetaNorm  Vector (not necessarily named) of booleans of the same dimension as the theta vector
#'   used to create the FIM (so not just the variable theta's).
#'   The value \code{TRUE} means the corresponding structural parameter should be normalized in the FIM,
#'   \code{FALSE} means not.
#'   Instead of a vector, may also provide a single boolean, where \code{TRUE} stands for a vector of
#'   all \code{TRUE}, and \code{FALSE} (default) for a vector of all \code{FALSE}.
#' @param omegaNorm  Vector (not necessarily named) of booleans of the same dimension as the omega matrix
#'   used to create the FIM (so not just the variable omega's).
#'   The value \code{TRUE} means the corresponding individual parameter should be normalized in the FIM,
#'   \code{FALSE} means not.
#'   Instead of a vector, may also provide a single boolean, where \code{TRUE} stands for a vector of
#'   all \code{TRUE}, and \code{FALSE} (default) for a vector of all \code{FALSE}.
#' @param sigmaNorm  Vector (not necessarily named) of booleans of the same dimension as the sigma matrix
#'   used to create the FIM (so not just the variable sigma's).
#'   The value \code{TRUE} means the corresponding residual parameter should be normalized in the FIM,
#'   \code{FALSE} means not.
#'   Instead of a vector, may also provide a single boolean, where \code{TRUE} stands for a vector of
#'   all \code{TRUE}, and \code{FALSE} (default) for a vector of all \code{FALSE}.
#' @param relChanges      If \code{FALSE} (default), the "directions" output contains absolute parameter changes.
#'   If \code{TRUE}, it contains parameter changes relative to the parameter values, as percentage, i.e., computed
#'   as 100 * (absolute change) / (parameter value). This option should not be combined with normalization, and
#'   a warning is printed if it is.
#' @param blockDiag       If \code{FALSE} (default), use the full FIM. If \code{TRUE}, use the block-diagonal FIM
#'   where mixed derivatives with respect to one structural and one random parameter are set to 0.
#' @param checkPosDefSymm If \code{TRUE}, it is checked that the FIM is positive semidefinite and symmetric.
#'   If \code{FALSE} (default), this is not checked.
#' @param checkPos        If \code{TRUE} (default), it is checked that the parameters to be normalized in the FIM
#'   are strictly positive. If \code{FALSE}, this is not checked.
#' @param ci              Confidence level. If \code{NULL} (default), the returned principal directions
#'   are vectors of length 1.
#'   If not \code{NULL}, they provide the parameter change corresponding to the given confidence level, see details.
#' @param nsubj          Number of subjects (default 1, should be > 0) used in the calculation of the curvature and
#'   the parameter change corresponding to the given confidence level.
#'
#' @return   \code{NULL} if there were errors. In this case an error message will be printed. 
#'   Otherwise, a list of seven elements, named "identifiable", "nDirections", "directions", "curvatures", "jump",
#'   "se", "rse", in that order.
#'   The first is \code{TRUE} in case the parameters are identifiable, and \code{FALSE} otherwise.
#'   The second contains the number of non-identifiable directions.
#'   The third contains the principal directions (eigenvectors of the FIM), in order of increasing curvature.
#'   This includes identifiable and non-identifiable directions.
#'   The direction vectors contain absolute or relative changes (depending on normalization and \code{relChanges})
#'   in parameter space.
#'   The fourth contains their curvatures (eigenvalues of the FIM), in the same order.
#'   The fifth contains the index of the curvature where the curvature increases by the largest factor.
#'   I.e. if the curvatures are c_1, ..., c_n, then it is the index i where c_{i+1}/c_i is the largest.
#'   If there are negative or zero curvatures, then it is set to the index of the last non-positive curvature.
#'   This could serve as a cutoff point between the unidentifiable (up to and including index "jump")
#'   and identifiable (higher than "jump") directions.
#'   The sixth contains the standard errors and the seventh the relative standard errors, as a percentage.
#'   If they cannot be computed (because the FIM cannot be inverted), then they are set to NA.
#'
#' @details The principal directions are the eigenvectors of the FIM, and their curvature is the corresponding eigenvalue.
#'   The eigenvector of a zero eigenvalue is a non-identifiable direction in parameter space,
#'   while positive eigenvalues correspond to identifiable directions.
#'   Small positive values would then be 'close to non-identifiable'.
#'   
#'   Normalization of the FIM can be used for parameters for which relative changes in value are the relevant quantity,
#'   rather than absolute changes, e.g., clearances, rates, volumes, and additive errors.
#'   Absolute changes in such parameters are not relevant because they depend on the chosen units.
#'   Relative changes do not suffer from this problem.
#'   It should typically also be used for parameters that already have a relative size,
#'   such as IIV variances (in case of multiplicative IIV effects, i.e., as exp(eta) factors),
#'   proportional errors, powers ("gamma"), Imax, etc.
#'   The reason is that these parameters may also be scaled by arbitrary multiplication factors
#'   (e.g., a proportional error can be modelled as a fraction or as a percentage).
#'   See \code{\link{normalizeFim}} for further details.
#'   Block-diagonalization is applied after normalization (not that it matters, the outcome is the same).
#'   
#'   The confidence interval and number of subjects are used to determine the amount by which the OFV can increase,
#'   assuming that it is chi-squared distributed with one degree of freedom.
#'   For each principal direction, this increase is achieved for a vector in this direction of a particular size.
#'   This vector is returned as principal vector.
#'   The only exception is when the curvature in a principal direction is zero or negative.
#'   Then the principal vector would have infinite length, and instead a vector of norm 1 is returned
#'   and a warning is given (if combined with \code{relChanges = TRUE} then the vector has norm 100%).
#'   The calculation assumes that each of the subjects contributes equally to the OFV.
#'   If this is not the case, then a separate FIM should be created for each subject, the sum of these FIMs
#'   should be provided as the parameter \code{fim}, and the number of subjects \code{nsubj} should be set to 1.
#'   
#'   Standard errors and relative standard errors are calculated from the diagonal of the inverse FIM
#'   (after applying normalization and block-diagonalization, if required), taking into account the number
#'   of subjects.
#'   
#' @export
#' 
#' @family calculation functions
#'
#' @author Martijn van Noort
fimIdent <- function(fim, curvature, thetaNorm = FALSE, omegaNorm = FALSE, sigmaNorm = FALSE, relChanges = FALSE,
                     blockDiag = FALSE, checkPosDefSymm = FALSE, checkPos = TRUE, ci = NULL, nsubj = 1) {
  # Check if relative direction is not combined with normalization:
  if (any(c(thetaNorm, omegaNorm, sigmaNorm)) && relChanges) {
    processWarns("FisherInfo::fimIdent: applying both normalization and relative parameter changes, which may give nonsensible results.\nContinuing.\n")
  }
  # Normalize FIM:
  nfim <- normalizeFim(fim = fim, thetaNorm = thetaNorm, omegaNorm = omegaNorm, sigmaNorm = sigmaNorm,
                       checkPos = checkPos)
  if(is.null(nfim)) return(NULL)  # There were errors
  # Optionally use the block-diagonal form:
  if (blockDiag) {
    thetaNams <- intersect(names(attr(nfim, thetaAttr)), colnames(nfim))
    randomNams <- intersect(c(colnames(attr(nfim, omegaAttr)), colnames(attr(nfim, sigmaAttr))), colnames(nfim))
    nfim[thetaNams, randomNams] <- nfim[randomNams, thetaNams] <- 0
  }
  # Compute eigenvectors and values, ordered from low to high.
  if (nsubj <= 0) {
    processErrors(paste0("FisherInfo::fimIdent: nsubj should be > 0.\nExiting.\n"))
    return(NULL)
  }
  es <- eigen(x = nfim, symmetric = TRUE)
  evals <- nsubj*rev(es$values)    # NB: this works better than multiplying the FIM one line above, due to num approx errors.
  evecs <- es$vectors
  evecs <- evecs[, rev(seq_len(ncol(evecs)))]
  row.names(evecs) <- row.names(nfim)
  # Check correctness
  if (checkPosDefSymm & any(evals < 0)) {
    processErrors(paste0("FisherInfo::fimIdent: FIM not positive semidefinite.\nExiting.\n"))
    return(NULL)
  }
  if (checkPosDefSymm & any(Im(evals) != 0)) {
    processErrors(paste0("FisherInfo::fimIdent: FIM not symmetric.\nExiting.\n"))
    return(NULL)
  }
  if (!is.null(ci)) {
    if (any(evals <= 0)) processWarns("FisherInfo::fimIdent: vectors will not be scaled for zero or negative curvatures.\nContinuing.\n")
    mult <- suppressWarnings(sqrt(qchisq(ci, 1)/evals))
    mult[is.na(mult) | is.infinite(mult)] <- 1
    evecs <- evecs %*% diag(mult, nrow = length(mult))
  }
  if (relChanges) {
    # Turn evecs into a set of relative changes:
    parms <- getParVecFim(fim)
    evecs <- 100 * evecs / parms   # This divides each row by the corresponding parameter value.
  }
  if (!is.null(ci) & relChanges) {
    # Use relative change of norm 100% in case curvature <= 0 and a CI is specified:
    evecnorms <- apply(evecs, MARGIN = 2, FUN = function(col) sqrt(col %*% col))
    evecs <- evecs %*% diag(ifelse(evals <= 0, 100/evecnorms, 1), nrow = length(evals))
  }
  # Compute index of largest jump in curvature:
  jump <- if(any(evals <= 0)) max(which(evals <= 0)) else which.max(diff(log(evals)))
  # Use jump <- which.max(diff(log(evals[evals > 0]))) + sum(evals <= 0) to find the largest jump among the positive ones.
  # Calculate SE's and RSE's:
  fiminv <- try(solve(nfim), silent = TRUE)
  if(inherits(fiminv, "try-error")) {
    se <- rep(NA, ncol(nfim))
    names(se) <- colnames(nfim)
    processWarns(paste0("FisherInfo::fimIdent: SE's and RSE's cannot be calculated because FIM cannot be inverted. I will set them to NA, and continue without them.\n"))
  } else {
    suppressWarnings(se <- sqrt(diag(fiminv)))
    if (any(is.na(se))) processWarns(paste0("FisherInfo::fimIdent: Some SE's and RSE's cannot be calculated because FIM is not positive definite. I will set them to NA, and continue without them.\n"))
  }
  pars <- getParVecFim(nfim)[names(se)]
  rse <- 100 * se / pars
  return(list(
    "identifiable" = all(evals >  curvature),
    "nDirections"  = sum(evals <= curvature),
    "directions"   = evecs,
    "curvatures"   = evals,
    "jump"         = jump,
    "se"           = se,
    "rse"          = rse
  ))
}


#------------------------------------------ simplifyFimIdent ------------------------------------------
#' Simplify FIM identifiability results
#'
#' Simplify the "directions" element of the output of \code{\link{fimIdent}}, by setting all elements
#' smaller than a threshold (in absolute value) to 0.
#' The resulting vectors are scaled to their original norms.
#'
#' @param fimId Named list of FIM identifiability indicators as produced by \code{\link{fimIdent}}.
#' @param tol  Threshold for setting elements to 0.
#'   Default value 0.001, i.e. all elements contributing less than 0.1% to the Euclidean vector norm are set to 0.
#'
#' @return   The list \code{fimId}, with a modified "directions" element.
#'   If the input list \code{fimId} does not contain this element, then it is returned without change.
#'
#' @export
#' 
#' @family result modifiers
#'
#' @author Martijn van Noort
simplifyFimIdent <- function(fimId, tol = 0.001) {
  if (!"directions" %in% names(fimId)) return(fimId)
  dirs <- fimId[["directions"]]
  nrms <- apply(dirs, MARGIN = 2, norm, type = "2")   # Euclidean norms of vectors in "directions"
  for (j in seq_along(nrms)) {
    dirs[, j] <- ifelse(abs(dirs[, j]) < tol * nrms[[j]], 0, dirs[, j])   # Set small elements of vector j to 0
    dirs[, j] <- dirs[, j] * nrms[[j]]/norm(dirs[, j], type = "2")        # Scale vector j to original norm.
  }
  fimId[["directions"]] <- dirs
  return(fimId)
}


#------------------------------------------ plotFim ------------------------------------------
#' Plot the FIM as a heatmap
#'
#' The given FIM is plotted as a heatmap, with one tile for each matrix element.
#' 
#' If on log scale, the heatmap uses separate color gradients for positive and negative values, if possible.
#' If this is not possible (it requires the package \code{ggnewscale}),
#' then negative values are indicated by a '-' symbol.
#' Optionally, the heatmap is labelled with the matrix values.
#'
#' @param fim   A matrix with row and column names containing the FIM,
#'   as produced by \code{\link{calcFimFromMatrix}} or \code{\link{calcFimFromModel}}.
#' @param vars  Vector of names of parameters to be included in the plots.
#'   If \code{NULL} (default) or missing, all parameters in the FIM are used.
#' @param log   If \code{TRUE} (default), then the FIM values are shown on the log scale.
#'   If \code{FALSE}, then not.
#' @param label If \code{TRUE}, then the heatmap is labelled with the FIM values.
#'   If \code{FALSE} (default), then not.
#'
#' @return   A plot.
#'
#' @export
#' 
#' @family plotting and printing
#'
#' @author Martijn van Noort
plotFim <- function(fim, vars = NULL, log = TRUE, label = FALSE) {
  if (is.null(vars)) vars <- row.names(fim)
  if (length(setdiff(vars, row.names(fim))) > 0) {
    processErrors("FisherInfo::plotFim: 'vars' contains parameters not present in 'fim'.\nExiting.\n")
    return(NULL)
  }
  # Create tall data set with FIM values versus the row and column indices i and j. Add signed values for log transforms:
  fim <- as.data.frame(fim)[vars, vars, drop = FALSE]
  fim[, "i"] <- factor(rownames(fim), levels = rownames(fim), ordered = TRUE)
  fim <- reshape2::melt(data = fim, id.vars = "i", variable.name = "j", value.name = "FIM")
  fim[, "j"] <- factor(as.character(fim[, "j"]), levels = rev(levels(fim[, "i"])), ordered = TRUE)
  fim[, "pos"] <- ifelse(fim[, "FIM"] > 0, fim[, "FIM"], NA)
  fim[, "neg"] <- ifelse(fim[, "FIM"] < 0, -1/fim[, "FIM"], NA)  # 1/ to reverse the order in the scale of the plot
  fim[, "sgn"] <- ifelse(fim[, "FIM"] < 0, "-", NA)
  # Package ggnewscale allows multiple scales of same type.
  # See https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/
  # and https://gist.github.com/eliocamp/eabafab2825779b88905954d84c82b32
  hasNewScale <- require(ggnewscale)
  clr <- rgb(0,0,0,alpha = 0)   # Fully transparent color for NA values. Use this instead of NA.value = NA, because that places tiles in the wrong locations (for some unknown reason)!
  pl <- 
    if (log & hasNewScale) {
      ggplot2::ggplot(fim) +
        geom_raster(aes(x = i, y = j, fill = pos), na.rm = TRUE) +
        scale_fill_gradient(element_blank(), low = "#132B43", high = "#56B1F7", trans = "log", na.value = clr,
                            guide = guide_colourbar(order = 1),
                            labels = function(brks) { formatC(x = brks, format = "g", digits = 3) }) +
        ggnewscale::new_scale("fill") +
        geom_raster(aes(x = i, y = j, fill = neg), na.rm = TRUE) +
        scale_fill_gradient(element_blank(), low = "violet", high = "violetred4", trans = "log", na.value = clr,
                            guide = guide_colourbar(order = 2),  # Always below the first one
                            labels = function(brks) { formatC(x = -1/brks, format = "g", digits = 3) })
    } else if (log) {
      pl1 <- ggplot2::ggplot(fim) +
        geom_raster(aes(x = i, y = j, fill = abs(FIM))) +
        scale_fill_gradient("Abs value", low = "#132B43", high = "#56B1F7", trans = "log", na.value = NA,
                            labels = function(brks) { formatC(x = brks, format = "g", digits = 3) })
      if (label) pl1 else pl1 + geom_text(aes(x = i, y = j, label = sgn), color = "white", na.rm = TRUE)
    } else {
      # Not on log scale.
      ggplot2::ggplot(fim) +
        geom_raster(aes(x = i, y = j, fill = FIM)) +
        scale_fill_gradient(element_blank(), low = "#132B43", high = "#56B1F7", na.value = NA,
                            labels = function(brks) { formatC(x = brks, format = "g", digits = 3) })
    }
  # Create title:
  titl <- NULL
  if (!label & !hasNewScale & log) titl <- c(titl, "- for negative values")
  if (hasNewScale & log & sum(fim[, "FIM"] == 0) > 0) titl <- c(titl, "empty for zero values")
  if (!is.null(titl)) titl <- paste0(" (", paste0(titl, collapse = ", "), ")")
  titl <- paste0("FIM", titl)
  pl <- pl + 
    scale_x_discrete("Parameter 1", position = "top") + ylab("Parameter 2") +
    ggtitle(titl) +
    ggplot2::theme_bw() + ggplot2::theme(legend.spacing = unit(0, "cm"), legend.margin = margin())
  if (label) pl <- pl + geom_text(aes(x = i, y = j, label = formatC(x = FIM, format = "g", digits = 3)), color = "white")
  return(pl)
}
