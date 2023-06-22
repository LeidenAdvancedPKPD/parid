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
#------------------------------------------ Exported functions ------------------------------------------
#--------------------------------------------------------------------------------------------------------


#------------------------------------------ calcSensitivityFromMatrix ------------------------------------------
#' Sensitivity indices, parameter relations and least identifiable parameters
#'
#' Calculates sensitivity indices, parameter relations and/or the least identifiable parameters for a given variational matrix.
#' This matrix may be parameter- and/or output-normalized.
#'
#' @param outputs    List of the indices and parameter relations to be determined.
#'   List of strings, where each element is one of "N", "S", "A", "R", "M", "L".
#'   They stand for null space dimension, skewing coefficient, skewing angle, parameter relations,
#'   minimal parameter relations and least identifiable parameters, respectively.
#' @param df         Data frame of variations (normalized or not), that is, with columns
#'   't', 'i', 'y' and 'dy_d(varsj)', where variables varsj are replaced by names.
#' @param vars    Vector of names of parameters to be included in sensitivity analysis.
#'   If \code{NULL} (default) or missing, all parameters in df are used.
#' @param Mthresh    Threshold for minimal parameter relations.
#'   Parameter relations smaller than \code{Mthresh/sqrt(k)} (before scaling, see \code{scaleM}) are returned,
#'   where k is the number of parameters.
#'   The factor sqrt(k) makes Mthresh independent of dimension.
#'   If \code{NULL} (default), only the most minimal relation is returned.
#' @param normalizeM \code{TRUE} (default) if the minimal parameter relations should be normalized, \code{FALSE} if not.
#'   This is not the same as parameter- or output-normalization of the sensitivity matrix.
#'   See the details for more information.
#' @param scaleM \code{TRUE} (default) if norms of the minimal parameter relations should be scaled to the range [0, 1],
#'   \code{FALSE} if not.
#'   If \code{TRUE}, the norms are divided by \code{sqrt(k)}, where k is the number of parameters.
#'   If \code{FALSE} they are unchanged.
#'
#' @return   A list with a named element for each element of 'outputs'. Elements for "N", "S", "A" are numeric values,
#'   while elements "R", "M" contain row-named matrices with the parameter relations defined by its columns,
#'   and "L" contains a named list with one element for each parameter, in increasing order of identifiability or distance,
#'   of the form "parameter name" = "distance of this parameter to the span of the others" in observation space.
#'   Distance (near) zero indicates non-identifiability of this parameter.
#'   Element M has an attribute "norm" specifying the norms of the minimal parameter relations, where near zero values
#'   correspond to non-identifiability, and large values to identifiability.
#'   The relations are provided in order of increasing norm (if there is more than one).
#'   It also has attributes "normalized" and "scaled" containing the values of the input parameters \code{normalizeM}
#'   and \code{scaleM}, respectively.
#'
#'   If \code{df} or \code{vars} is not defined correctly, then the return value is \code{NULL}
#'   and an error message is shown.
#'   If there is a single parameter in \code{vars}, then elements S, A and L are \code{NA}.
#'
#' @details  
#' The sensitivity calculations are based on the sensitivity matrix df = dy/dp, i.e.,
#' the derivative of the model output y with respect to the parameters p.
#' Non-trivial vectors in the null space of this matrix indicate directions in which the parameter can be changed
#' (infinitesimally) without changing the model output.
#' A null space of dimension > 0 therefore indicates potential non-identifiability, while a trivial null space
#' implies identifiability.
#'
#' The skewing coefficient (S) measures the extent to which the parameters are mapped parallel to
#' each other in the observations y.
#' The coefficient can take any non-negative value.
#' A small coefficient indicates that the model is (close to) non-identifiable.
#' For a single parameter the coefficient is by definition equal to 1.
#' A disadvantage is that the value of this indicator depends on the number of parameters.
#'
#' The skewing angle (A) is a dimensionless version (i.e. independent of the number of parameters)
#' of the skewing coefficient that measures the angle (in radians) between the images of
#' parameter vectors in the space of observations y.
#' A small angle indicates that the model is (close to) non-identifiable.
#' For a single parameter the angle is by definition equal to 1.
#'
#' Minimal parameter relations (M) are identified as the parameter vector(s) with the smallest image
#' (under the sensitivity matrix or a normalized version of it) in the space of observations.
#' The image is calculated under a normalized version of the sensitivity matrix.
#' This can be used in case of non-identifiability or identifiability (but is most useful near non-identifiability).
#' With normalization (i.e., \code{normalizeM = TRUE}), this indicator finds the combination of parameters for which
#' the directions of their images are closest to cancelling each other out, without taking into account the sizes of
#' the image vectors.
#' Without normalization, these sizes are taken into account, and the minimal parameter relation would typically
#' point more towards the parameter that least influences the output.
#' Both variants map the same parameter vectors to 0, if any, so they agree on the distinction between
#' non-identifiability and identifiability, but may differ regarding what is nearly non-identifiable.
#'
#' Finally, the least identifiable parameter(s) (L) can be identified as the parameters whose image in the space of
#' observations are nearest to the span of the images of the other parameters.
#' The images are calculated under a normalized version of the sensitivity matrix.
#' This can be used in case of non-identifiability or identifiability (but is most useful near non-identifiability).
#'
#' @export
#' 
#' @family calculation functions
#'
#' @author Martijn van Noort
calcSensitivityFromMatrix <- function(outputs, df, vars = NULL, Mthresh = NULL, normalizeM = TRUE, scaleM = TRUE) {
  isVarMat <- isValidVariations(df)
  if (!is.null(isVarMat)) {
    processErrors(paste0("Sensitivity::calcSensitivityFromMatrix: 'df' is not a variational matrix.\n", paste0(paste0(unlist(isVarMat), "\n"), collapse = ""), "Exiting.\n"))
    return(NULL)
  }
  df <- df[, grep("^d2y", value = TRUE, invert = TRUE, names(df))]
  if (is.null(vars)) vars <- gsub("^dy_d", "", grep("^dy_d", names(df), value = TRUE))
  colms <- paste0("dy_d", vars)
  nonExist <- setdiff(colms, names(df))
  if(length(nonExist) > 0) {
    processErrors("Sensitivity::calcSensitivityFromMatrix: 'vars' contains parameters not present in variational matrix.\nExiting.\n")
    return(NULL)
  }

  # Obtain matrix dy/dp (or a normalized version):
  mx        <- as.matrix(df[, colms, drop = FALSE])
  npars     <- ncol(mx)
  out       <- list()

  # "N":
  nullspace <- if ("N" %in% outputs | "R" %in% outputs) {
    res <- MASS::Null(t(mx))
    row.names(res) <- vars
    res
  } else {
    NA
  }
  if ("N" %in% outputs) out <- c(out, list("N" = ncol(nullspace)))

  # "S":
  skew <- if ("S" %in% outputs | "A" %in% outputs) {
    if (npars <= 1) {
      1
    } else {
      stretch   <- prod(apply(mx, 2, function(col) sqrt(sum(col**2))))   # Amount of 'stretching' by dy/dp.
      sqmx      <- t(mx) %*% mx
      # volUnit   <- if (ncol(MASS::Null(sqmx)) > 0) 0 else sqrt(abs(det(sqmx))) # Volume under dy/dp of unit ball. Determinant does not detect non-invertibility well, so use MASS::Null for that
      volUnit   <- sqrt(abs(det(sqmx))) # Volume under dy/dp of unit ball.
      if (stretch > 0) volUnit / stretch else 0   # If stretch = 0, it means that dy/dp_i = 0 for some parameter p_i, and hence unidentifiable
    }
  } else {
    NA
  }
  if ("S" %in% outputs) out <- c(out, list("S" = skew))

  # "A":
  if ("A" %in% outputs) {
    angle <- if (npars <= 1) {
      1
    } else {
      skew**(1/(npars-1))
    }
    out <- c(out, list("A" = angle))
  }

  # "R":
  if ("R" %in% outputs) {
    out       <- c(out, list("R" = nullspace))
  }

  norms <- if ("M" %in% outputs | "L" %in% outputs) {
    norms1 <- apply(X = mx, MARGIN = 2, FUN = function(col) sqrt(col %*% col))  # column norms of dY/dp
    ifelse(norms1 == 0, 1, norms1)
  } else {
    NA
  }

  normmx <- if ("M" %in% outputs | "L" %in% outputs) {
    scale(mx, center = FALSE, scale = norms)     # normmx = dY/dp with columns normalized
  } else {
    NA
  }

  # "M":
  if ("M" %in% outputs) {
    mx1                          <- if(normalizeM) normmx else mx
    eigensys                     <- eigen(t(mx1) %*% mx1)  # Eigensystem of mx^t mx or normmx^t normmx
    # Find indices of selected eigenvectors:
    selInds                      <- rev(if(is.null(Mthresh)) which.min(eigensys$values) else which(eigensys$values < Mthresh**2/ncol(mx1)))
    eigenvecs                    <- eigensys$vectors[, selInds, drop = FALSE]  # Selected eigenvectors, in order of increasing eigenvalues
    minrels                      <- diag(1/norms, nrow = length(norms)) %*% eigenvecs
    minrelNrm                    <- apply(X = minrels, MARGIN = 2, FUN = function(col) sqrt(col %*% col))  # column norms
    minrels                      <- scale(minrels, center = FALSE, scale = minrelNrm)                      # normalize columns
    row.names(minrels)           <- vars
    attr(minrels, "norm")        <- sqrt(pmax(eigensys$values[selInds], 0))/(if(scaleM) sqrt(ncol(mx1)) else 1)  # Min norms = sqrt(Eigenvalues) of selected Eigenvectors
    attr(minrels, "normalized")  <- normalizeM
    attr(minrels, "scaled")      <- scaleM
    attr(minrels,"scaled:scale") <- NULL
    out                          <- c(out, list("M" = minrels))
  }

  # "L":
  if ("L" %in% outputs & npars > 1) {
    remainders <- unlist(lapply(1:ncol(normmx), function(i) {
      # remainders[[i]] is the norm of the component of normmx[, i] orthogonal to the span of the columns of normmx[, -i].
      ui <- normmx[, i, drop = FALSE]    # i-th column from normalized dY/dp
      Vi <- normmx[, -i, drop = FALSE]   # normalized dY/dp with i-th column dropped
      wi <- try(ui - Vi %*% solve(t(Vi) %*% Vi, t(Vi) %*% ui), silent = TRUE)  # difference between ui and its projection onto span(Vi).
      if (inherits(wi, "try-error")) {
        # Presumably the inversion failed because columns of Vi are dependent. Create independent subset:
        s <- c()   # This will be the indices of independent vectors in Vi
        for (j in 1:ncol(Vi)) {
          Vij <- Vi[, j, drop = FALSE]
          for (k in s) {
            Vij <- Vij - c(t(Vij) %*% Vi[, k, drop = FALSE]) * Vi[, k, drop = FALSE]  # Vij = component of Vi[, j] orthogonal to span of columns in Vi[, s]
            # Vij[abs(Vij) < 1e-10] <- 0  # Not necessary because of use of all.equal below, which tests near equality.
          }
          if (isTRUE(all.equal(Vij, 0*Vij))) {
            # Orthogonal component is 0, so Vij does not add to span of columns of Vi
            Vi[, j] <- 0
          } else {
            # Orthogonal component is not 0, so Vij adds to span of columns of Vi
            Vi[, j] <- Vij / c(sqrt(t(Vij) %*% Vij))  # This makes columns of Vi orthonormal, which is required for GramSchmidt like procedure used above.
            s <- c(s, j)
          }
        }
        Vi <- Vi[, s, drop = FALSE]  # Restrict to independent columns s.
        wi <- ui - Vi %*% solve(t(Vi) %*% Vi, t(Vi) %*% ui)
      }
      return(sqrt(t(wi) %*% wi))
    }))
    names(remainders) <- vars                              # names of parameters
    ord       <- order(remainders)                         # ordering by increasing remainder
    out       <- c(out, list("L" = remainders[ord]))
  } else if ("L" %in% outputs) out <- c(out, list("L" = NA))

  # Output:
  return(out)
}

#------------------------------------------ calcSensitivityFromModel ------------------------------------------
#' Sensitivity indices, parameter relations and least identifiable parameters
#'
#' Calculates sensitivity indices, parameter relations and/or the least identifiable parameters for a given model.
#' It first computes the variational matrix for the model and then applies \code{\link{calcSensitivityFromMatrix}}.
#'
#' @param outputs    List of the indices and parameter relations to be determined.
#'   List of strings, where each element is one of "N", "S", "A", "R", "M", "L".
#'   They stand for null space dimension, skewing coefficient, skewing angle, parameter relations,
#'   minimal parameter relations and least identifiable parameters, respectively.
#' @param model      Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
#'   This function should return the numeric vector dy/dt.
#' @param parms      Model parameter values, as named numeric vector.
#' @param init       Function(p) of parameters specifying the initial state of the model.
#' @param outputPred Function(y, p) of state variables and parameters specifying the model prediction, i.e., the
#'   output without residual.
#'   This function should return a numeric vector.
#' @param times      Numeric vector of times where variations are to be evaluated.
#' @param vars       Vector of names of variables q to be included in sensitivity analysis.
#'   Should be a subset of \code{names(parms)}.
#'   By default equal to \code{names(parms)}.
#'   Variations d(outputPred)/dq are calculated only for these parameters.
#'   If specified as \code{NULL}, it is set to \code{names(parms)}.
#' @param parNorm    Vector (not necessarily named) of booleans of the same dimension as \code{parms}.
#'   The value \code{TRUE} means the corresponding individual parameter should be normalized in the FIM,
#'   \code{FALSE} means not.
#'   If the vector is named, the names should match \code{vars}.
#'   Instead of a vector, may also provide a single boolean, where \code{TRUE} stands for a vector of
#'   all \code{TRUE}, and \code{FALSE} (default) for a vector of all \code{FALSE}.
#' @param outNorm    \code{TRUE} if variations should be output-normalized, \code{FALSE} (default) if not.
#' @param symbolic   \code{TRUE} (default) if derivatives are to be computed symbolically, \code{FALSE} if numerically.
#'   See \code{\link{calcVariations}} for notes on using symbolic derivation.
#' @param chkModel   \code{TRUE} (default) if it has to be checked whether model components (\code{model}, \code{p}, \code{init}, \code{output}) 
#'   are formatted correctly for symbolic manipulation, \code{FALSE} if not.
#' @param Mthresh    Threshold for minimal parameter relations.
#'   Parameter relations smaller than \code{Mthresh/sqrt(k)} (before scaling, see \code{scaleM}) are returned,
#'   where k is the number of parameters.
#'   The factor sqrt(k) makes Mthresh independent of dimension.
#'   If \code{NULL}, only the most minimal relation is returned.
#' @param normalizeM \code{TRUE} (default) if the minimal parameter relations should be normalized, \code{FALSE} if not.
#'   This is not the same as parameter- or output-normalization of the sensitivity matrix.
#'   See \code{\link{calcSensitivityFromMatrix}} for more information.
#' @param scaleM \code{TRUE} (default) if norms of the minimal parameter relations should be scaled to the range [0, 1],
#'   \code{FALSE} if not.
#'   If \code{TRUE}, the norms are divided by \code{sqrt(k)}, where k is the number of parameters.
#'   If \code{FALSE} they are unchanged.
#' @param ...        Named arguments to be passed to \code{\link[deSolve]{lsoda}}.
#'   Can be used for example to pass events or error tolerances.
#'
#' @return   A list with a named element for each element of 'outputs'. Elements for "N", "S", "A" are numeric values,
#'   while elements "R", "M" contain row-named matrices with the parameter relations defined by its columns,
#'   and "L" contains a named list with one element for each parameter, in increasing order of identifiability or distance,
#'   of the form "parameter name" = "distance of this parameter to the span of the others" in observation space.
#'   Distance (near) zero indicates non-identifiability of this parameter.
#'   Element M has an attribute "norm" specifying the norms of the minimal parameter relations, where near zero values
#'   correspond to non-identifiability, and large values to identifiability.
#'   The relations are provided in order of increasing norm (if there is more than one).
#'   It also has attributes "normalized" and "scaled" containing the values of the input parameters \code{normalizeM}
#'   and \code{scaleM}, respectively.
#'
#'   If the variational matrix could not be computed or \code{vars} is not defined correctly,
#'   then the return value is \code{NULL} and an error message is shown.
#'   If there is a single parameter in \code{vars}, then elements S, A and L are \code{NA}.
#'
#' @seealso \code{\link{calcSensitivityFromMatrix}} for details on the calculation of the outputs.
#'
#' @note Output-normalization divides the vector dy_i/dp_j by y_i, and is recommended if the output is
#'   considered to be a multiplicative rather than an additive variable.
#'   Parameter-normalization multiplies the vector dy_i/dp_j by p_j, and is recommended if the input is
#'   considered to be a multiplicative rather than an additive variable.
#'   (A multiplicative variable is a variable where differences between values are considered on the log scale.
#'   So if e.g. values 0.001 and 0.01 are considered different because they differ by a factor of 10, then the
#'   variable is multiplicative.
#'   If they are considered approximately equal because they are both close to 0,
#'   then the variable is additive.
#'   Most PK inputs and outputs are multiplicative, likewise for many PD ones.)
#'   
#' @note Output-normalization is equivalent to considering the sensitivity matrix dlog(y)/dp rather than dy/dp.
#'   Parameter-normalization is equivalent to considering the sensitivity matrix dy/dlog(p) rather than dy/dp.
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
calcSensitivityFromModel <- function(outputs, model, parms, init, outputPred, times, vars = names(parms),
                                     parNorm = FALSE, outNorm = FALSE, symbolic = TRUE, chkModel= TRUE,
                                     Mthresh = NULL, normalizeM = TRUE, scaleM = TRUE, ...) {
  # This simply calculates the variations and then calls the companion method that works on the variation matrix in tall format.
  if (is.null(vars)) vars <- names(parms)
  df <- calcFirstVariations(model = model, parms = parms, init = init, outputPred = outputPred, times = times,
                            varp = vars, symbolic = symbolic, chkModel = chkModel, ...)
  if (is.null(df)) return(NULL)
  df <- normalizeVariations(df = df, parNorm = parNorm, outNorm = outNorm)
  if (is.null(df)) return(NULL)
  return(calcSensitivityFromMatrix(outputs = outputs, df = df, vars = vars, Mthresh = Mthresh,
                                   normalizeM = normalizeM, scaleM = scaleM))
}

#------------------------------------------ simplifySensitivities ------------------------------------------
#' Simplify sensitivity results
#'
#' Simplify the "R", "M" and/or "L" element of the output of \code{\link{calcSensitivityFromMatrix}} or
#' \code{\link{calcSensitivityFromModel}}, by setting all elements smaller than a threshold (in absolute value) to 0.
#' The results are renormalized.
#'
#' @param sens Named list of sensitivity indicators as produced by \code{\link{calcSensitivityFromMatrix}} or
#'   \code{\link{calcSensitivityFromModel}}.
#' @param elt  List with elements "R", "M" and/or "L" for simplification of R, M and/or L.
#'   \code{NULL} can be used as shorthand for all.
#' @param tol  Threshold for setting elements to 0.
#'   Default value 0.001, i.e. all elements contributing less than 0.1% to the Euclidean vector norm are set to 0.
#'
#' @return   The list \code{sens}, with a modified "R", "M" and/or "L" element.
#'   If the input list \code{sens} does not contain these elements, then it is returned without change.
#'
#' @export
#' 
#' @family result modifiers
#'
#' @author Martijn van Noort
simplifySensitivities <- function(sens, elt = NULL, tol = 0.001) {
  allCases <- c("R", "M", "L")
  if (is.null(elt)) elt <- allCases
  for (el in elt) {
    if (el %in% names(sens) & el %in% allCases) {
      if (el == "M") dist <- attr(sens[[el]], "norm")  # preserve "norm" attribute
      sens[[el]] <- ifelse(abs(sens[[el]]) < tol, 0, sens[[el]])
      if (el %in% c("R", "M")) {
        # Normalize columns of matrix:
        norms <- apply(sens[[el]], MARGIN = 2, FUN = function(col) sqrt(col %*% col))
        sens[[el]] <- scale(sens[[el]], center = FALSE, scale = norms)
        attr(sens[[el]], "scaled:scale") <- NULL
      }
      if (el == "M") attr(sens[[el]], "norm") <- dist  # preserve "norm" attribute
    }
  }
  return(sens)
}

#------------------------------------------ plotSensitivities ------------------------------------------
#' Plot sensitivity results
#'
#' Generates plots for one or more of the "R", "M" and "L" elements of the output of
#' \code{\link{calcSensitivityFromMatrix}}, \code{\link{calcSensitivityFromModel}} or \code{\link{simplifySensitivities}}.
#' For "R" or "M", the plot shows the parameter relations, while for "L" it shows the ordering of the parameters.
#'
#' @param sens Named list of sensitivity indicators as produced by \code{\link{calcSensitivityFromMatrix}},
#'   \code{\link{calcSensitivityFromModel}} or \code{\link{simplifySensitivities}}.
#' @param elt  List with elements "R", "M" and/or "L" indicating which plot(s) should be made.
#'   \code{NULL} can be used as shorthand for all.
#'
#' @return A named list of plots with one element for each element in elt that is present in sens.
#'   In case of error, the return value is \code{NULL}.
#'
#' @export
#' 
#' @family plotting and printing
#'
#' @author Martijn van Noort
plotSensitivities <- function(sens, elt = NULL) {
  allCases <- c("R", "M", "L")
  if (is.null(elt)) elt <- allCases else elt <- intersect(elt, allCases)
  elt <- intersect(elt, names(sens))
  pls <- lapply(elt, function(el) {
    if (el %in% c("R", "M")) {
      df <- as.data.frame(sens[[el]])
      colNames <- colnames(df)
      nVars <- length(colNames)
      if(nVars == 0) return(NULL)
      df[, "par"] <- factor(row.names(df), levels = row.names(df), ordered = TRUE)
      nPars <- nrow(df)
      df <- reshape2::melt(df, id.vars = "par")
      out <- ggplot2::ggplot(df) + ggplot2::theme_bw() + ylab("Weight") +
        scale_x_continuous("Parameter", breaks = 1:nPars, labels = levels(df[, "par"]))
      if(nPars > 1) out <- out + geom_vline(data = data.frame(par = (0:nPars)+0.5), mapping = aes(xintercept = par), color = "lightgrey")
      if(el == "R") {
        return(out +
                 geom_col(aes(x = as.numeric(par), y = value, fill = variable), position = "dodge", show.legend = FALSE) +
                 ggtitle(paste0("Parameter relation", if(nVars > 1) "s (each one a different color)" else ""))
        )
      }
      if(el == "M") {
        norms <- format(attr(sens[[el]], "norm"), digits = 2)
        names(norms) <- colNames
        cols <- grDevices::colorRampPalette(c("#000060", "#0000FF"))(nVars)   # From dark to lightblue
        names(cols) <- colNames
        return(out +
                 geom_col(aes(x = as.numeric(par), y = value, fill = variable), position = "dodge") +
                 scale_fill_manual("Norm", values = cols, breaks = colNames, labels = norms) +
                 ggtitle(paste0("Minimal parameter relation", if(nVars > 1) "s" else ""))
        )
      }
    }
    if (el == "L") {
      if (is.na(sens[[el]][[1]])) return(NULL)
      df <- data.frame(par = factor(names(sens[[el]]), levels = names(sens[[el]]), ordered = TRUE), value = sens[[el]])
      cols <- grDevices::colorRampPalette(c("#000060", "#0000FF"))(nrow(df))   # From dark to lightblue
      return(ggplot2::ggplot(df) + geom_col(aes(x = par, y = value, fill = par), show.legend = FALSE) + ggplot2::theme_bw() +
               scale_fill_manual("", values = cols, breaks = levels(df[, "par"])) +
               ggtitle("Parameters in order from least to most identifiable") + xlab("Parameter") + ylab("Distance"))
    }
  })
  names(pls) <- elt
  return(pls[!is.null(pls)])
}
