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


#------------------------------------------ calcAliasingScoresFromMatrix ------------------------------------------
#' Normalized sensitivities and aliasing scores
#'
#' Calculate the normalized sensitivities and aliasing scores for a given variational matrix.
#' This matrix may be parameter- and/or output-normalized.
#' Parameter-normalization will not change the outcome (as the sensitivities are normalized anyway), but output-normalization may do so.
#'
#' @param df   Data frame of variations in tall format (normalized or not), that is, with columns
#'   't', 'i', 'y' and 'dy_d(varsj)', where variables varsj are replaced by names.
#' @param vars Vector of names of parameters to be included in aliasing analysis. If \code{NULL} or missing, all parameters in df are used.
#'
#' @details    
#' Aliasing scores are calculated as follows.
#' The sensitivity of the model to parameter p_j is s_j(t,i) = dy_i/d_p_j(t), where y_i is the i-th output variable and t is time.
#' The normalized sensitivity is sn_j(t,i) = s_j(t,i) / max_{t,i}|s_j(t,i)|, and takes values between -1 and 1.
#' The aliasing score of of parameters p_j and p_k is defined as a(j,k) = 100% * (1-max_{t,i}(| |sn_j(t,i)| - |sn_k(t,i)| |)),
#' i.e. it measures the similarity between the profiles of their normalized sensitivities.
#' In particular a(j,j) = 100%.
#' For more details see Augustin et al (2019) and Braakman (2019).
#'
#' @references Augustin, Braakman, Paxson - A workflow to detect non-identifiability in parameter estimation
#'   using SimBiology, Github (2019).
#'   Sietse Braakman - A general workflow for parameter estimation to help establish confidence in model
#'   predictions, Rosa Webinar 2019-12-01, slides 39-42.
#'
#' @return     A list with two elements, named "normSens" and "aliasingScore".
#'   The first contains a matrix with as named columns the time "t", output indicator index "i",
#'   and normalized sensitivities of the parameters selected by \code{vars}.
#'   The second contains a matrix with the aliasing scores, one for each pair of parameters
#'   from the selection in \code{vars}.
#'   Each score is a percentage between 0 and 100, where 100 indicates similarity between the parameters
#'   and unidentifiability, and 0 means the parameters are not similar, indicating identifiability.
#'   The matrix is symmetric and the diagonal elements are 100.
#'   In case of errors, the return value is \code{NULL}.
#'
#' @export
#'
#' @author Martijn van Noort
calcAliasingScoresFromMatrix <- function(df, vars = NULL) {
  isVarMat <- isValidVariations(df)
  if (!is.null(isVarMat)) {
    processErrors(paste0("Aliasing::calcAliasingScoresFromMatrix: 'df' is not a variational matrix.\n", paste0(paste0(unlist(isVarMat), "\n"), collapse = ""), "\nExiting.\n"))
    return(NULL)
  }
  df <- df[, grep("^d2y", value = TRUE, invert = TRUE, names(df))]  # Remove second order derivs as we do not need them.
  if (is.null(vars)) vars <- gsub("^dy_d", "", grep("^dy_d", names(df), value = TRUE))
  colms <- paste0("dy_d", vars)
  nonExist <- setdiff(colms, names(df))
  if(length(nonExist) > 0) {
    processErrors("Aliasing::calcAliasingScoresFromMatrix: 'vars' contains parameters not present in variational matrix.\nExiting.\n")
    return(NULL)
  }

  # Obtain matrix dy/dp (or a normalized version):
  mx        <- as.matrix(df[, colms, drop = FALSE])
  ind       <- df[, "i"]   # output index
  npars     <- ncol(mx)
  if (npars < 1) return(NULL)  # Nothing to do

  # Calculate normalized sensitivities:
  nmx <- t(plyr::aaply(.data = mx, .margins = 2, .drop = FALSE, .fun = function(col) {
    mxcol <- tapply(X = abs(col), INDEX = ind, FUN = max)  # Calculate maximum value per output index.
    mxcol <- ifelse(mxcol > 0, mxcol, 1)                   # Set 0 values to 1 (no normalization in this case)
    out <- col / mxcol[as.character(ind)]                  # Divide each element by the appropriate maximum value
    names(out) <- NULL                                     # Otherwise names are set to "i" and get warning "Duplicate names not supported".
    out
  }))  # Transpose because aaply permutes the dimensions

  # Calculate pairwise aliasing scores:
  aliasing <- matrix(data = 100, nrow = npars, ncol = npars, dimnames = list(i = vars, j = vars))
  for (i in seq(1, npars, by = 1)) {
    if (npars > i) {
      for (j in seq(i+1, npars, by = 1)) {
        aliasing[[i,j]] <- 100*max(0, 1-max(abs(abs(nmx[, i]) - abs(nmx[, j]))))
      }
    }
    if (i > 1) for (j in seq(1, i-1, by = 1)) aliasing[[i,j]] <- aliasing[[j,i]]
  }

  # Output:
  return(list(normSens = cbind(df[, c("t", "i")], nmx), aliasingScore = aliasing))
}



#------------------------------------------ calcAliasingScoresFromModel ------------------------------------------
#' Normalized sensitivities and aliasing scores
#'
#' Calculate the normalized sensitivities and aliasing scores for a given model.
#' It first computes the variational matrix for the model and then applies \code{\link{calcAliasingScoresFromMatrix}}.
#'
#' @param model      Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
#'   This function should return the numeric vector dy/dt.
#' @param parms      Model parameter values, as named numeric vector.
#' @param init       Function(p) of parameters specifying the initial state of the model.
#' @param outputPred Function(y, p) of state variables and parameters specifying the model prediction, i.e., the
#'   output without residual.
#'   This function should return a numeric vector.
#' @param times      Numeric vector of times where variations are to be evaluated.
#' @param vars       Vector of names of variables q to be included in aliasing analysis.
#'   Should be a subset of \code{names(parms)}.
#'   By default equal to \code{names(parms)}.
#'   Variations d(output)/dq are calculated only for these parameters.
#'   If specified as \code{NULL}, it is set to \code{names(parms)}.
#' @param outNorm    \code{TRUE} if variations should be output-normalized, \code{FALSE} (default) if not.
#' @param symbolic   \code{TRUE} (default) if derivatives are to be computed symbolically, \code{FALSE} if numerically.
#'   See \code{\link{calcVariations}} for notes on using symbolic derivation.
#' @param chkModel   \code{TRUE} (default) if it has to be checked whether model components (\code{model}, \code{p}, \code{init}, \code{output}) 
#'   are formatted correctly for symbolic manipulation, \code{FALSE} if not.
#' @param ...        Named arguments to be passed to \code{\link[deSolve]{lsoda}}.
#'   Can be used for example to pass events or error tolerances.
#'
#' @return A list with two elements, named "normSens" and "aliasingScore".
#'   The first contains a matrix with as named columns the time "t", output indicator index "i",
#'   and normalized sensitivities of the parameters selected by \code{vars}.
#'   The second contains a matrix with the aliasing scores, one for each pair of parameters
#'   from the selection in \code{vars}.
#'   Each score is a percentage between 0 and 100, where 100 indicates similarity between the parameters
#'   and unidentifiability, and 0 means the parameters are not similar, indicating identifiability.
#'   The matrix is symmetric and the diagonal elements are 100.
#'   In case of errors, the return value is \code{NULL}.
#'
#' @seealso \code{\link{calcAliasingScoresFromMatrix}} for details on the calculation of the scores.
#'
#' @note Output-normalization divides the vector dy_i/dp_j by y_i, and is recommended if the output is
#'       considered to be a multiplicative rather than an additive variable.
#'       (A multiplicative variable is a variable where differences between values are considered on the log scale.
#'       So if e.g. values 0.001 and 0.01 are considered different because they differ by a factor of 10, then the
#'       variable is multiplicative.
#'       If they are considered approximately equal because they are both close to 0,
#'       then the variable is additive.
#'       E.g. most PK outputs are multiplicative.)
#'       
#' @note Output-normalization is equivalent to considering the sensitivity matrix dlog(y)/dp rather than dy/dp.
#' 
#' @note Parameter-normalization would not change the outcome (as the sensitivities are normalized anyway),
#'       and hence there is no input parameter \code{parNorm}.
#'       
#' @note It is recommended to use symbolic derivation (\code{symbolic==TRUE}) rather than numerical derivation,
#'   as it is more accurate, especially for derivatives that are formally zero.
#'   See \code{\link{calcVariations}} for details.
#'
#' @export
#'
#' @author Martijn van Noort
calcAliasingScoresFromModel <- function(model, parms, init, outputPred, times, vars = names(parms),
                                        outNorm = FALSE, symbolic = TRUE, chkModel = TRUE, ...) {
  # This simply calculates the variations and then calls the companion method that works on the variation matrix in tall format.
  if (is.null(vars)) vars <- names(parms)
  df <- calcFirstVariations(model = model, parms = parms, init = init, outputPred = outputPred, times = times, varp = vars,
                            symbolic = symbolic, chkModel = chkModel, ...)
  if (is.null(df)) return(NULL)
  df <- normalizeVariations(df = df, parNorm = FALSE, outNorm = outNorm)
  if (is.null(df)) return(NULL)
  return(calcAliasingScoresFromMatrix(df, vars))
}


#------------------------------------------ plotAliasing ------------------------------------------
#' Plot of normalized sensitivities and aliasing scores
#'
#' Generates several plots for normalized sensitivities and/or aliasing scores, as computed by
#' \code{\link{calcAliasingScoresFromMatrix}} or \code{\link{calcAliasingScoresFromModel}}.
#'
#' @param aliasing Named list of aliasing score output as produced by \code{\link{calcAliasingScoresFromMatrix}}
#'   or \code{\link{calcAliasingScoresFromModel}}.
#' @param elt      List with elements "S", "P", "A" and/or "T" indicating which plot(s) should be made.
#'   \code{NULL} can be used as shorthand for all.
#' @param vars     Vector of names of parameters to be included in the plots.
#'   If \code{NULL} (default) or missing, all parameters in aliasing are used.
#'
#' @return A named list of plots with the same names as in elt.
#'   Element "S" is a time evolution plot of the normalized sensitivities.
#'   Element "P" is a time evolution plot of the absolute value of the normalized sensitivities.
#'   Element "A" is a heatmap of the aliasing scores.
#'   Element "T" is like "A", but also showing the aliasing scores as text.
#'
#' @export
#'
#' @author Martijn van Noort
plotAliasing <- function(aliasing, elt = NULL, vars = NULL) {
  if (is.null(vars)) vars <- colnames(aliasing[["aliasingScore"]])
  colms <- paste0("dy_d", vars)
  nonExist <- setdiff(vars, colnames(aliasing[["aliasingScore"]]))
  if(length(nonExist) > 0) {
    processErrors("Aliasing::plotAliasing: 'vars' contains parameters not present in 'aliasing'.\nExiting.\n")
    return(NULL)
  }
  allCases <- c("S", "P", "A", "T")
  if (is.null(elt)) elt <- allCases else elt <- intersect(elt, allCases)
  pls <- lapply(elt, function(el) {
    if (el %in% c("S", "P")) {
      # Plot time evolutions of Normalized sensitivities
      df <- as.data.frame(aliasing[["normSens"]])[, c("t", "i", colms)]
      df <- reshape2::melt(df, id.vars = c("t", "i"))
      df[, "i"] <- as.character(df[, "i"])
      df[, "variable"] <- factor(as.character(df[, "variable"]), levels = colms, ordered = TRUE)
      if (el == "P") df[, "value"] <- abs(df[, "value"])
      txt <- if (el == "P") "Absolute n" else "N"
      pl <- ggplot2::ggplot(df) +
        geom_line(aes(x = t, y = value, group = interaction(variable, i), color = variable, linetype = i)) +
        ggplot2::theme_bw() + xlab("Time") + ylab(paste0(txt, "ormalized sensitivity")) +
        ggtitle(paste0(txt, "ormalized sensitivities over time")) +
        scale_color_discrete("")
      if (length(unique(df[, "i"])) > 1) pl <- pl + labs(linetype = "Output nr") else pl <- pl + guides(linetype = "none")
      return(pl)
    } else if (el %in% c("A", "T")) {
      # Plot heatmap or table of Aliasing scores
      df <- as.data.frame(aliasing[["aliasingScore"]])[vars, vars, drop = FALSE]
      df[, "i"] <- factor(rownames(df), levels = rownames(df), ordered = TRUE)
      df <- reshape2::melt(data = df, id.vars = "i", variable.name = "j", value.name = "Score")
      df[, "j"] <- factor(as.character(df[, "j"]), levels = rev(levels(df[, "i"])), ordered = TRUE)
      df[df[, "i"] == df[, "j"], "Score"] <- NA
      df[, "ScoreTxt"] <- paste0(round(df[, "Score"]), "%")
      pl <- ggplot2::ggplot(df) + geom_raster(aes(x = i, y = j, fill = Score)) +
        ggplot2::theme_bw() + scale_x_discrete("Parameter 1", position = "top") + ylab("Parameter 2") +
        scale_fill_continuous("Score (%)", limits = c(0, 100), breaks = seq(0, 100, by = 25), na.value="white") +
        ggtitle("Aliasing score (high = similar)")
      if (el == "T") {
        pl <- pl + geom_text(aes(x = i, y = j, label = ScoreTxt), color = "white")
      }
      return(pl)
    }
  })
  names(pls) <- elt
  return(pls)
}

