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


#------------------------------------------ calcVariationsNum ------------------------------------------
#' Calculate variations numerically
#'
#' Calculates the variational matrix of the given model, that is, the first and optionally second order derivatives of
#' the model outputs with respect to the model parameters, at given times.
#'
#' @param model    Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
#'   This function should return the numeric vector dy/dt.
#' @param p        Function(theta, eta) of population and individual parameters specifying the model parameters.
#'   This function should return a named numeric vector.
#' @param init     Function(p) of parameters specifying the initial state (at time 0) of the model.
#' @param output   Function(y, p, eps) of state variables, parameters and residual errors specifying the model outputs.
#'   This function should return a numeric vector.
#' @param times    Numeric vector of times where variations are to be evaluated.
#'   Should all be >= 0, but do not need to include the dose time 0.
#' @param theta    Population parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param eta      Individual parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param eps      Residual parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param vartheta Vector of names of population parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(theta)}, or \code{NULL} for none.
#' @param vareta   Vector of names of individual parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(eta)}, or \code{NULL} for none.
#' @param vareps   Vector of names of residual parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(eps)}, or \code{NULL} for none.
#' @param secOrd   \code{TRUE} if second order derivatives have to be computed, \code{FALSE} if not.
#' @param ...    Named arguments to be passed to \code{\link[deSolve]{lsoda}}.
#'   Can be used for example to pass events or error tolerances.
#'
#' @return A data frame with columns 't' for time, 'i' for the index of the output element, 'y' for the outputs,
#'   'dy_d<v1>' for the first derivatives and (if \code{secOrd == TRUE}) 'd2y_d<v1>_<v2>' for the second derivatives.
#'   In the column names, variables v1 and v2 are replaced by the names in \code{vartheta}, \code{vareta} and \code{vareps}.
#'   The function displays an error and returns \code{NULL} if integration of the model fails, or if \code{vartheta}
#'   contains elements not in \code{names(theta)}, and likewise for \code{eta} and \code{eps}.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
calcVariationsNum <- function(model, p, init, output, times, theta, eta, eps,
                              vartheta, vareta, vareps, secOrd, ...) {
  # Create the complements of vartheta, vareta and vareps, and the corresponding parameter vectors:
  for (nam in c("theta", "eta", "eps")) {
    assign(paste0("var", nam, "c"), setdiff(names(get(nam)), get(paste0("var", nam))))   # Defines varthetac, varetac, varepsc as complements of vartheta, vareta, vareps.
    assign(paste0("val", nam), get(nam)[get(paste0("var", nam))])                        # Defines value vectors valtheta  = theta[vartheta], etc.
    assign(paste0("val", nam, "c"), get(nam)[get(paste0("var", nam, "c"))])              # Defines value vectors valthetac = theta[varthetac], etc.
  }
  # Now valtheta, valeta, valeps contain the parameters for which derivatives are to be computed;
  # and valthetac, valetac, valepsc contain the ones that stay fixed.
  
  # Create the ODE and output function to be used in lsoda:
  fu <- function(t, y, p, eps) {
    c(list(model(t, y, p)), output(y, p, eps))
  }
  parms <- p(theta, eta)
  nOutputs <- length(output(init(parms), parms, eps))
  xTimes <- sort(unique(c(0, times)))
  # Create the function that relates the input parameters to the output, as a data frame:
  varfu <- function(q) {
    # q contains valtheta, valeta,valeps
    theta1  <- c(q[vartheta], valthetac)
    eta1    <- c(q[vareta],   valetac)
    eps1    <- c(q[vareps],   valepsc)
    parms   <- p(theta1, eta1)
    out <- try(deSolve::lsoda(y = init(parms), times = xTimes, func = fu, parms = parms, eps = eps1, ...))
    if(inherits(out, "try-error")) return(NULL)
    return(out[out[, "time"] %in% times, (ncol(out)-nOutputs+1):ncol(out), drop = FALSE])
    # Keep the final nOutputs columns corresponding to the output, as a matrix.
    
    # return(as.numeric(out[, (ncol(out)-nOutputs+1):ncol(out), drop = FALSE]))
    # # Keep the final nOutputs columns corresponding to the output, and list them as single vector. (As single vector is not needed.)
    # # Can go back to matrix form with 'matrix(vector, ncol = nOutputs)'.
  }
  # Calculate the derivatives:
  q <- c(valtheta, valeta, valeps)
  derivs <- if(is.null(q) || length(q) == 0) list(varfu(q), NULL) else if(secOrd) numDeriv::genD(varfu, q)[c(3, 1)] else list(varfu(q), numDeriv::jacobian(varfu, q))
  # This produces a list with two elements:
  # The first is a matrix with one column per output element, and one row per time point.
  # The second is a matrix of derivatives, with one column per derivative, and the output element columns below each other as rows (or NULL if there were no elements in q).
  # Now create output. First a data frame with times and simulated point values:
  out <- data.frame(derivs[[1]])
  names(out) <- 1:nOutputs    #paste0("y", 1:nOutputs)
  out <- cbind(data.frame(t = times), out)
  out <- reshape2::melt(out, id.vars = "t", variable.name = "i", value.name = "y")
  out[, "i"] <- as.numeric(as.character(out[, "i"]))   # One column each for time, output index, and output value.
  # Then the derivatives:
  if(is.null(derivs[[2]])) return(out)    # No derivatives to add
  derivs <- data.frame(derivs[[2]])
  nams <- paste0("dy_d", names(q))   # names of first derivs
  if (secOrd) {
    # add names of second derivs
    nams2 <- matrix(paste0("d2y_d", names(q), "_", rep(names(q), each = length(q))), nrow = length(q))
    nams <- c(nams, nams2[upper.tri(nams2, diag = TRUE)])
  }
  names(derivs) <- nams
  derivs[, "t"] <- rep(times, times = nOutputs)
  derivs[, "i"] <- rep(1:nOutputs, each = length(times))
  # Return combined results:
  return(merge(out, derivs, all = TRUE, by = c("t", "i")))
}

#------------------------------------------ copyEnv ------------------------------------------
#' Copy an environment
#'
#' Creates a shallow copy of the given environment.
#' Changes to this copy will not affect the original environment.
#' In this sense, \code{newEnv <- copyEnv(env)} is different from \code{newEnv <- env}.
#' Designed to create a copy of the rules environment of \code{\link[Deriv]{Deriv}}.
#'
#' @param env    Environment to copy.
#'
#' @return A copy of this environment.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
copyEnv <- function(env) {
  as.environment(as.list(env, all.names=TRUE))
}


#------------------------------------------ createLowerTriSecondDerivCodes ------------------------------------------
#' Create second derivative codes
#'
#' Creates second derivative codes from first derivative ones in lower triangular form, that is, codes are generated
#' for derivatives (i,j) where i <= j.
#' 
#' @param firstDerivCodes    Vector of strings containing codes of first derivatives.
#' @param sep                String defining the separator character(s) between the two derivatives in the output.
#'   Default is ".".
#'
#' @return Vector of strings containing second derivative codes, of the form \code{<deriv1><sep><deriv2>}, where
#'   \code{<deriv1>} and \code{<deriv2>} are elements of \code{firstDerivCodes}.
#'   Returns \code{NULL} if there are none. This happens when \code{firstDerivCodes} is empty or NULL.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
createLowerTriSecondDerivCodes <- function(firstDerivCodes, sep = ".") {
    df <- expand.grid(j = seq_along(firstDerivCodes), i = seq_along(firstDerivCodes))
    df <- df[df[, "i"] <= df[, "j"], ]
    if (nrow(df) == 0) NULL else paste0(firstDerivCodes[df[, "i"]], sep, firstDerivCodes[df[, "j"]])
}
  

#------------------------------------------ createDerivGreps ------------------------------------------
#' Create grep strings for derivatives.
#'
#' Creates grep strings for selecting the given derivatives from the output of \code{\link[Deriv]{Deriv}}.
#' To this end, the codes provided in \code{derivCodes} are translated to match the names of elements of this output.
#' The grep strings take into account that these names can be different, depending on simplifications of the output.
#' That is, output vectors equal to 0 are sometimes reduced to single elements, and corresponding parts of their
#' name are suppressed.
#' Also, derivatives may be wholly or partially suppressed.
#'
#' @param derivCodes    Vector of strings containing the derivative codes.
#'   Elements have the form \code{0} for the 0-th derivative, \code{1<derivSep><deriv>}
#'   for the 1-st derivatives, and \code{2<derivSep><deriv><derivSep><deriv>} for the second
#'   derivatives. Here \code{<deriv>} codes a derivative as \code{<type><typeSep><variable>},
#'   where \code{<type>} is a derivative type, such as "theta" or "eta", and \code{<variable>}
#'   is a variable (from the list of variables of the given type), such as "TVCL".
#' @param compNames     Vector of strings specifying the names of the output components of the 0-th derivative,
#'                      or \code{NULL} if they have no names.
#' @param varName       String specifying the variable name to be used in the output.
#' @param derivSep      String specifying the character(s) separating the various parts of an element of \code{derivCodes}.
#' @param typeSep       String specifying the character(s) separating the type and variable in elements of \code{derivCodes}.
#' @param outnameSep    String specifying the character(s) separating the various parts of the derivative specification
#'   in the names of the output.
#'                      
#' @return The grep strings, as a named vector of strings of the same length as derivCodes.
#'   The names of its elements are \code{<varName><outnameSep><deriv1><outnameSep><deriv2>}, where the derivatives are
#'   included as needed.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
createDerivGreps <- function(derivCodes, compNames, varName, derivSep = "\\.", typeSep = "_", outnameSep = "_") {
  splitCodes <- strsplit(derivCodes, split = derivSep)
  splitCodesNames <- unlist(lapply(splitCodes, function(splitCode) paste0(c(varName, gsub(paste0("^.*", typeSep), "", splitCode[-1])), collapse = "_")))
  typeOpt <- FALSE   # Is type part of name optional? Can set to TRUE if it is.
  splitCodes <- unlist(lapply(splitCodes, function(splitCode) {
    n <- length((splitCode))
    splitCode[-1] <- gsub(paste0("(.*)", typeSep, "(.*)"), paste0("(\\\\.(\\1_)", if(typeOpt) "?" else "", "\\2"), splitCode[-1])
    compCode <- if(is.null(compNames)) "[[:digit:]]*$" else paste0("(|\\.(", paste0(compNames, collapse = "|"), "))$") # Either a number (possibly NULL) or one of the given names preceded by a "."
    closingCode <- paste0(rep(")?", n-1), collapse = "")
    paste0("^", paste0(splitCode, collapse = ""), closingCode, compCode)
  }))
  names(splitCodes) <- splitCodesNames
  return(splitCodes)
}


#------------------------------------------ transformExpr ------------------------------------------
#' Transform an expression by replacing variables by functions or vice versa
#' 
#' Transforms an expression (e.g. a function body) by replacing state variables by state functions or vice versa.
#' The state variables are assumed to be of the form 'y[[<state>]]', where state is a unique identifier
#' of the component, of the form '<Component>_<Deriv1>_<Deriv2>', where the component may be a name or a number,
#' and the derivatives are optional.
#' 
#' @param expr     The expression to transform.
#' @param pList    List of parameters to be used as function arguments. Each element is a language object.
#' 
#' @return The transformed expression.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
transformExpr <- function(expr, pList) {
  lpList <- length(pList)
  if (is.atomic(expr) || is.name(expr)) {
    return(expr)
  } else if (is.call(expr)) {
    # Check if expr equals y[[<state>]]:
    isStateElt <- length(expr) >= 3 &&
      identical(expr[[1]], quote(`[[`)) &&
      identical(expr[[2]], quote(y))
    # Check if expr equals FU_<state_deriv> (this is a reserved function name used to code state variables):
    isStateFunction <- length(expr) == 1+lpList &&
      grepl("^FU_", as.character(expr[[1]]))
    for (i in 1:lpList) isStateFunction <- isStateFunction && identical(pList[[i]], expr[[i+1]])
    if (isStateElt) {
      # Change expr to function if isStateElt:
      expr <- as.call(c(as.name(paste0("FU_", expr[[3]])), pList))
    } else if (isStateFunction) {
      # Change expr to y[[]] if isStateFunction:
      expr <- as.call(c(quote(`[[`), quote(y), gsub("^FU_", "", expr[[1]])))
    } else expr <- as.call(unlist(lapply(X = expr, FUN = transformExpr, pList = pList))) # Call transformExpr recursively for other expressions
    return(expr)
  } else if (is.pairlist(expr)) {
    return(c("[", unlist(lapply(X = expr, FUN = transformExpr, pList = pList)), "]")) # Call transformExpr recursively for pairlists
  } else {
    # User supplied incorrect input
    stop("Don't know how to handle type ", typeof(expr), 
         call. = FALSE)
  }
}

    
#------------------------------------------ createDerivsSymb ------------------------------------------
#' Generate derivatives of model components symbolically
#'
#' Calculates first and optionally second derivatives of the model, initialization, output and parameter functions
#' with respect to relevant parameters, in a symbolic way.
#'
#' @param model    Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
#'   This function should return the numeric vector dy/dt.
#' @param p        Function(theta, eta) of population and individual parameters specifying the model parameters.
#'   This function should return a named numeric vector.
#' @param init     Function(p) of parameters specifying the initial state (at time 0) of the model.
#' @param output   Function(y, p, eps) of state variables, parameters and residual errors specifying the model outputs.
#'   This function should return a numeric vector.
#'   
#' @param theta    Population parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param eta      Individual parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param eps      Residual parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param vartheta Vector of names of population parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(theta)}, or \code{NULL} for none.
#' @param vareta   Vector of names of individual parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(eta)}, or \code{NULL} for none.
#' @param vareps   Vector of names of residual parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(eps)}, or \code{NULL} for none.
#' @param secOrd   \code{TRUE} if second order derivatives have to be computed, \code{FALSE} if not.
#'
#' @return A list of four functions, containing the variational versions of \code{model}, \code{p}, \code{init} and
#'   \code{output}. These functions take same arguments as their non-variational counterparts, but with a different
#'   meaning for the argument \code{y}, namely this contains the state variables and their derivatives wrt parameters.
#'   The output of \code{p} is a named list of named numeric vectors, with one element for each derivative.
#'   The names of the elements have the format \code{p(_<deriv>(_<deriv>))}, where \code{<deriv>} codes the derivative and
#'   is optional, as indicated by the brackets.
#'   Each element has the same length as the output of the non-variational version of \code{p}, and its components have
#'   the same names. 
#'   The output of \code{init} and \code{model} is a named numeric vector of values of the state variables and their
#'   derivatives.
#'   The names have the format \code{y(_<deriv>(_<deriv>)).<state>}, where \code{<deriv>} codes the derivative and
#'   is optional, and \code{<state>} codes the state.
#'   The output of \code{output} is a named numeric vector of values of the output variables and their derivatives.
#'   The names have the format \code{y(_<deriv>(_<deriv>)).<state>}, where \code{<deriv>} codes the derivative and
#'   is optional, and \code{<state>} codes the state or (if the state has no name) its index.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
createDerivsSymb <- function(model, p, init, output, theta, eta, eps, vartheta, vareta, vareps, secOrd) {
  # ---
  # Preliminaries.
  # Create outputs of model components (this is used to determine their format):
  pval <- p(theta, eta)
  ival <- init(pval)
  mval <- model(0, ival, pval)
  oval <- output(ival, pval, eps)
  # Create flags that determine whether model components are referred by name (TRUE) or by index (FALSE):
  pHasnames <- length(names(pval)) == length(pval)
  iHasnames <- length(names(ival)) == length(ival)
  mHasnames <- length(names(mval)) == length(mval)
  oHasnames <- length(names(oval)) == length(oval)
  # Create names of model component output elements. If there are no names, use index as name:
  namps <- if(length(names(pval)) != length(pval)) seq_along(pval) else names(pval)
  namys <- if(length(names(ival)) != length(ival)) seq_along(ival) else names(ival)   # Used as y values (just like in lsoda)
  namms <- if(length(names(mval)) != length(mval)) seq_along(mval) else names(mval)   # Only used to parse derivative expression
  namos <- if(length(names(oval)) != length(oval)) seq_along(oval) else names(oval)
  # Determine which orders of derivatives are to be computed:
  nderiv <- if(secOrd) 0:2 else 0:1
  # Create named input parameters (theta, eta, eps):
  if(length(vartheta) > 0) names(vartheta) <- rep("theta", length(vartheta))
  if(length(vareta) > 0)   names(vareta)   <- rep("eta", length(vareta))
  if(length(vareps) > 0)   names(vareps)   <- rep("eps", length(vareps))
  varthetaeta     <- c(vartheta, vareta)
  
  # ---
  # Create dp = function of (theta, eta) that contains the derivatives of p(theta, eta).
  if(is.null(varthetaeta)) {
    # No derivatives (wrt theta or eta):
    dp <- eval(substitute(function(theta, eta) list('p' = f(theta, eta), 'p' = NULL, 'p' = NULL)[seq_along(nd)],
                          list(f = p, nd = nderiv)))
  } else {
    # There are derivatives (wrt theta or eta).
    # Inputs for derivative of p: dpCodes, dpNeutrVec, dpare:
    dpCodes <- {
      dpPars <- paste(names(varthetaeta), varthetaeta, sep = "_")   # Derivatives wrt these pars.
      dpPars2 <- createLowerTriSecondDerivCodes(dpPars)
      createDerivGreps(derivCodes = c("0", paste0("1.", dpPars), if (secOrd) paste0("2.", dpPars2) else NULL),
                       compNames = if(pHasnames) namps else NULL,
                       varName = "p")
    }
    dpNeutrVec <- rep(1, length(namps))
    parStr <- paste0(deparse(body(p)), collapse = "")
    dpare  <- Deriv::Deriv(parse(text = parStr), x = varthetaeta, nderiv = nderiv)
    
    # Derivative of p:
    dp <- eval(substitute(function(theta, eta) {
      evalExpr <- unlist(eval(expr))
      # Arrange in one vec per derivative:
      lapply(dCodes, function(dCode) {
        out1 <- evalExpr[grep(dCode, names(evalExpr))] * neutrVec  # Multiply by neutrVec in case vector length was reduced to 1
        names(out1) <- nams
        return(out1)
      })
    }, list(nams = namps, dCodes = dpCodes, neutrVec = dpNeutrVec, expr = dpare)))
  }
  
  # Select parameter names (from "p") that have nonzero derivatives. In model, init and output, only their derivatives are needed:
  nampsNonzeroDeriv <- namps[unlist(lapply(namps, function(namp) {
    any(unlist(lapply(dp(theta, eta)[-1], function(dpComp) {
      dpComp[[namp]] != 0
    })))
  }))]
  names(nampsNonzeroDeriv) <- rep("p", length(nampsNonzeroDeriv))
  if(length(nampsNonzeroDeriv) == 0) nampsNonzeroDeriv <- NULL
  
  # ---
  # Create dmodel = function of (t, y, p) that contains the derivatives of the ODE's model(t, y, p).
  # Initialize mydrule which contains additional rules for derivation of state variables ("FU_"...) and pList which
  # contains relevant elements of the parameter vector (p), as language objects.
  # Both are used in dmodel, doutput:
  mydrule <- copyEnv(Deriv::drule)
  pList <- list()
  if(is.null(nampsNonzeroDeriv)) {
    # No derivatives (wrt p):
    dmodel <- eval(substitute(function(t, y, p) {
      out <- f(t, y, p)
      names(out) <- nams
      return(out)
    }, list(nams = namys, f = model)))
  } else {
    # There are derivatives (wrt p).
    # Inputs for derivative of model (pvec, pList, FU_<x> dummy functions, dmodele, dmodelCodes, dmodelNeutrVec):
    pvec  <- if(pHasnames) paste0("p[[\"", nampsNonzeroDeriv, "\"]]") else paste0("p[[", nampsNonzeroDeriv, "]]")
    pList <- lapply(pvec, function(txt) parse(text = txt)[[1]])  # List with elements of parameter vector, as language objects
    pvec  <- paste0(pvec, collapse = ", ")
    # Create state functions for 0th, 1st and 2nd derivs, and rules for the derivation of 0th and 1st derivs:
    fargs <- rep(0, length(nampsNonzeroDeriv))
    names(fargs) <- nampsNonzeroDeriv
    for (namy in namys) {
      # Functions for 0th derivs:
      fnamy <- paste0("FU_", namy)  # Use prefix "FU_" to code state variables as functions.
      assign(fnamy, as.function(as.list(c(fargs, NaN))))
      # rules for taking 1st deriv:
      rules <- lapply(nampsNonzeroDeriv, function(namp) parse(text = paste0(fnamy, "_", namp, "(", pvec, ")"))[[1]])
      names(rules) <- nampsNonzeroDeriv
      mydrule[[fnamy]] <- rules
      for (i in seq_along(nampsNonzeroDeriv)) {  # i <- 1
        # Functions for 1st derivs:
        dfnamy <- paste0(fnamy, "_", nampsNonzeroDeriv[[i]])
        assign(dfnamy, as.function(as.list(c(fargs, NA))))
        # rules for taking 2nd deriv:
        drules <- lapply(seq_along(nampsNonzeroDeriv), function(j) {
          namp1 <- nampsNonzeroDeriv[[min(i, j)]]
          namp2 <- nampsNonzeroDeriv[[max(i, j)]]
          as.symbol(paste0(fnamy, "_", namp1, "_", namp2))
          parse(text = paste0(fnamy, "_", namp1, "_", namp2, "(", pvec, ")"))[[1]]
        })
        names(drules) <- nampsNonzeroDeriv
        mydrule[[dfnamy]] <- drules
        for (j in seq_along(nampsNonzeroDeriv)) {  # j <- 1
          # Functions for 2nd derivs:
          if (j >= i) {
            dfnamy <- paste0(fnamy, "_", nampsNonzeroDeriv[[i]], "_", nampsNonzeroDeriv[[j]])
            assign(dfnamy, as.function(as.list(c(fargs, NA))))
          }
        }
      }
    }
    dmodele <- transformExpr(body(model), pList)  # model with state variables replaced by functions
    dmodele <- Deriv::Deriv(f = dmodele, x = nampsNonzeroDeriv, nderiv = nderiv, drule = mydrule)  # Take derivatives
    dmodele <- transformExpr(dmodele, pList)     # Transform back to state variables
    dmodelCodes <- {
      dmodelPars <- paste(names(nampsNonzeroDeriv), nampsNonzeroDeriv, sep = "_")   # Derivatives wrt model pars.
      dmodelPars2 <- createLowerTriSecondDerivCodes(dmodelPars)
      createDerivGreps(derivCodes = c("0", paste0("1.", dmodelPars), if(secOrd) paste0("2.", dmodelPars2) else NULL),
                       compNames = if(mHasnames) namms else NULL,
                       varName = "y")
    }  # grep strings for derivatives
    dmodelNeutrVec <- rep(1, length(namys))
    
    # Derivative of model:
    dmodel <- eval(substitute(function(t, y, p) {
      evalExpr <- unlist(eval(expr))
      out <- unlist(lapply(dCodes, function(dCode) {
        out1 <- evalExpr[grep(dCode, names(evalExpr))] * neutrVec  # Multiply by neutrVec in case vector length was reduced to 1
        names(out1) <- nams
        return(out1)
      }))
      names(out) <- gsub("^y(.*)\\.(.*)$", "\\2\\1", names(out))  # To form <Component>_<Deriv1>_<Deriv2>, where derivs are optional
      return(out)
    }, list(nams = namys, dCodes = dmodelCodes, neutrVec = dmodelNeutrVec, expr = dmodele)))
  }
  
  # ---
  # Create dinit = function of (p) that contains the derivatives of the initial value init(p).
  if(is.null(nampsNonzeroDeriv)) {
    # No derivatives (wrt p):
    dinit <- eval(substitute(function(p) {
      out <- f(p)
      names(out) <- nams
      return(out)
    }, list(nams = namys, f = init)))
  } else {
    # There are derivatives (wrt p).
    # Inputs for derivative of init (dinite):
    dinite <- Deriv::Deriv(f = body(init), x = nampsNonzeroDeriv, nderiv = nderiv)  # Take derivatives
    dinitCodes <- {
      dinitPars <- paste(names(nampsNonzeroDeriv), nampsNonzeroDeriv, sep = "_")   # Derivatives wrt model pars.
      dinitPars2 <- createLowerTriSecondDerivCodes(dinitPars)
      createDerivGreps(derivCodes = c("0", paste0("1.", dinitPars), if(secOrd) paste0("2.", dinitPars2) else NULL),
                       compNames = if(iHasnames) namys else NULL,
                       varName = "y")
    }  # grep strings for derivatives
    
    # Derivative of init:
    dinit <- eval(substitute(function(p) {
      evalExpr <- unlist(eval(expr))
      out <- unlist(lapply(dCodes, function(dCode) {
        out1 <- evalExpr[grep(dCode, names(evalExpr))] * neutrVec  # Multiply by neutrVec in case vector length was reduced to 1
        names(out1) <- nams
        return(out1)
      }))
      names(out) <- gsub("^y(.*)\\.(.*)$", "\\2\\1", names(out))  # To form <Component>_<Deriv1>_<Deriv2>, where derivs are optional
      return(out)
    }, list(nams = namys, dCodes = dinitCodes, neutrVec = dmodelNeutrVec, expr = dinite)))
  }
  
  # ---
  # Create doutput = function of (y, p, eps) that contains the derivatives of the output function output(y, p, eps).
  if(is.null(nampsNonzeroDeriv) && is.null(vareps)) {
    # No derivatives (wrt p or eps):
    doutput <- eval(substitute(function(y, p, eps) {
      out <- f(y, p, eps)
      names(out) <- nams
      return(out)
    }, list(nams = namos, f = output)))
  } else {
    # There are derivatives (wrt p or eps).
    if(is.null(nampsNonzeroDeriv)) {
      # Derivatives wrt eps but not wrt p.
      # This needs to be treated differently because mydrule does not contain any rules for derivatives of y components.
      # And it can be treated differently, because only derivatives with respect to eps are needed:
      doutpute <- Deriv::Deriv(f = body(output), x = vareps, nderiv = nderiv)
    } else {
      # There are derivatives (wrt p or eps).
      # Inputs for derivative of output (doutpute, doutputCodes, doutputNeutrVec):
      doutpute <- transformExpr(body(output), pList)  # output with state variables replaced by functions
      doutpute <- Deriv::Deriv(f = doutpute, x = c(nampsNonzeroDeriv, vareps), nderiv = nderiv, drule = mydrule)  # Take derivatives
      doutpute <- transformExpr(doutpute, pList)     # Transform back to state variables
    }
    doutputCodes <- {
      doutputPars <- paste(names(c(nampsNonzeroDeriv, vareps)), c(nampsNonzeroDeriv, vareps), sep = "_")   # Derivatives wrt model pars and residuals.
      doutputPars2 <- createLowerTriSecondDerivCodes(doutputPars)
      createDerivGreps(derivCodes = c("0", paste0("1.", doutputPars), if(secOrd) paste0("2.", doutputPars2) else NULL),
                       compNames = if(oHasnames) namos else NULL,
                       varName = "y")
    }  # grep strings for derivatives
    doutputNeutrVec <- rep(1, length(namos))
    
    # Derivative of output:
    doutput <- eval(substitute(function(y, p, eps) {
      evalExpr <- unlist(eval(expr))
      # exprNams <- if(length(nams) <= 1) names(evalExpr) else stringr::str_sub(names(evalExpr), 1, -2)  # Take substring to remove index of output element; 
      out <- unlist(lapply(dCodes, function(dCode) {
        #      out1 <- evalExpr[grep(dCode, exprNams)] * neutrVec  # Multiply by neutrVec in case vector length was reduced to 1
        out1 <- evalExpr[grep(dCode, names(evalExpr))] * neutrVec  # Multiply by neutrVec in case vector length was reduced to 1
        names(out1) <- nams
        return(out1)
      }))
      names(out) <- gsub("^y(.*)\\.(.*)$", "\\2\\1", names(out))  # To form <Component>_<Deriv1>_<Deriv2>, where derivs are optional
      return(out)
    }, list(nams = namos, dCodes = doutputCodes, neutrVec = doutputNeutrVec, expr = doutpute)))
  }
  
  # Return the functions defining the variational equations:
  return(list("model" = dmodel, "p" = dp, "init" = dinit, "output" = doutput))
}


#------------------------------------------ extractDf ------------------------------------------
#' Extract columns from a data frame
#'
#' Extracts the specified columns from a data frame.
#' The difference with \code{[]} is that missing columns are replaced by zeroes.
#'
#' @param df    data frame
#' @param cols  vector of strings specifying columns to be extracted
#' 
#' @return A data frame containing the required columns
#' 
#' @author Martijn van Noort
#' 
#' @noRd 
extractDf <- function(df, cols) {
  out <- as.data.frame(matrix(0, nrow = nrow(df), ncol = length(cols)))
  names(out) <- cols
  for (col in cols) if(col %in% names(df)) out[, col] <- df[, col]
  return(out)
}


#------------------------------------------ calcVarSymb ------------------------------------------
#' Calculate output in case of symbolic derivatives
#'
#' Calculates the variations in case of symbolic derivatives, for a given output and parameter vector, both in variational form,
#' that is, including derivatives.
#' For the output vector these are derivatives with respect to the derived parameters (\code{p}).
#' For the parameter vector these are derivatives with respect to the base parameters (\code{theta}, \code{eta}, \code{eps}),
#' denoted here as \code{v}.
#' This function combines these to obtain derivatives of the output with respect to the base parameters.
#' This is an applicaton of the chain rule for derivatives: dy_i/dv = sum_j(dy_i/dp_j dp_j/dv) and
#' d^2y_i/(dv dw) = sum_{j,k}(d^2y_i/(dp_j dp_k) dp_j/dv dp_k/dv) + sum_j(dy_i/dp_j d^2p_j/(dv dw)).
#'
#' @param doutput    data frame containing values of derivatives of order 0, 1 and optionally 2 of the output
#'   vector, each row in the format of the output of the function \code{doutput}, see \code{\link{createDerivsSymb}},
#'   and each row corresponding to a time point.
#'   The column names should reflect the component and derivative, as in the output of the function \code{doutput}.
#' @param dpout     named numeric list containing values of derivatives of order 0, 1 and optionally 2 of the parameter
#'   vector, in the format of the output of the function \code{dp}, see \code{\link{createDerivsSymb}}.
#' @param vartheta Vector of names of population parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(theta)}, or \code{NULL} for none.
#' @param vareta   Vector of names of individual parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(eta)}, or \code{NULL} for none.
#' @param vareps   Vector of names of residual parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(eps)}, or \code{NULL} for none.
#' @param secOrd   \code{TRUE} if second order derivatives have to be computed, \code{FALSE} if not.
#'
#' @return A vector with derivatives of order 0, 1, 2 of the state variables with respect to base parameters
#'   (theta, eta and eps). The elements are named <nam>.<i>, where i is the sequence number of the component of y, 
#'   and nam is of the form 'y', 'dy_d<v1>' or 'd2y_d<v1>_<v2>'.
#'   
#' @author Martijn van Noort
#' 
#' @noRd 
calcVarSymb <- function(doutput, dpout, vartheta, vareta, vareps, secOrd) {
  namos <- unique(unlist(lapply(names(doutput), function(nam) strsplit(nam, "_")[[1]][[1]])))  # Component names (or indices) of output variables, in non-variational context
  namps <- names(dpout[[1]])  # Names of derived parameters
  namesHessian <- {
    df <- expand.grid(i = seq_along(namps), j = seq_along(namps))
    df[, "name"] <- paste0(namps[pmin(df[, "i"], df[, "j"])], "_", namps[pmax(df[, "i"], df[, "j"])])
    matrix(df[, "name"], nrow = length(namps), dimnames = list(namps, namps))
  }
  derivCases <- {
    dv  <- c(vartheta, vareta, vareps)
    dv2 <- createLowerTriSecondDerivCodes(dv, sep = "_")
    rep(c("", dv, if(secOrd) dv2 else NULL), each = length(namos))
  }  # vector of strings containing the derivatives to be taken ("" for no derivative, "var" for first
  #    derivative in var, "var1_var2" for second derivative in var1 and var2).
  out <- mapply(derivCase = derivCases, namo = namos, SIMPLIFY = FALSE, function(derivCase, namo) {
    derivs  <- strsplit(derivCase, "_")[[1]]  # Derivatives are taken with respect to these base parameters (may be empty)
    derivs  <- derivs[!derivs %in% vareps]    # Restrict to theta and eta derivatives
    derivNr <- length(derivs)                 # Number of theta and eta derivatives
    if (derivNr == 0) return(doutput[, paste0(namo, ifelse(nchar(derivCase) == 0, "", "_"), derivCase)])
    if (derivNr == 1) {
      namp <- paste0("p_", derivs)
      dpout1 <- dpout[[namp]]
      doutput1 <- extractDf(doutput, paste0(namo, "_", unlist(lapply(namps, function(namp) gsub(derivs, namp, derivCase)))))  # The construction "unlist(...)" is needed to cover also cases with residual derivatives
      return(unlist(plyr::alply(doutput1, .margins = 1, function(row) as.numeric(row) %*% dpout1)))
    }
    if (derivNr == 2) {
      namp1 <- paste0("p_", derivs)     # Names of first derivs
      namp2 <- paste0("p_", derivCase)  # Name of second deriv
      dpout1 <- dpout[namp1]            # List of two vectors of first derivs of p wrt namp1
      dpout2 <- dpout[[namp2]]          # Vector of second derivs of p wrt namp2
      doutput1 <- extractDf(doutput, paste0(namo, "_", namps))   # Since derivNr = 2, there are no derivatives wrt residual, so no need to account for them.
      doutput2 <- lapply(1:nrow(doutput), function(i) apply(namesHessian, MARGIN = c(1,2), function(elt) as.numeric(extractDf(doutput[i, ], paste0(namo, "_", elt)))))
      out1 <- unlist(lapply(doutput2, function(doutput2e) as.numeric(dpout1[[1]] %*% doutput2e %*% dpout1[[2]])))
      out2 <- unlist(plyr::alply(doutput1, .margins = 1, function(row) as.numeric(row) %*% dpout2))
      return(out1+out2)
    }
  })
  out <- data.frame(out)
  nDerivs <- as.character(unlist(lapply(strsplit(derivCases, "_"), length)))
  nDerivCodes <- c("0" = "y", "1" = "dy_d", "2" = "d2y_d")
  names(out) <- paste0(nDerivCodes[nDerivs], derivCases, ".", seq_along(namos))
  return(out)
}


#------------------------------------------ analyzeBody ------------------------------------------
#' Analyze structure of a function body
#'
#' Dissects a given function body (technically, any expression) to determine if forbidden functions \code{with},
#' \code{return} or \code{FU_...} occur anywhere, and to determine how the state \code{y} and parameter
#' \code{p} are indexed.
#' Functions \code{FU_...} are used internally to represent the state variables and therefore may not be used
#' for other purposes.
#' It does these checks by going recursively through the expression tree.
#' This function assumes that \code{y} and \code{p} are the formal arguments representing the state and
#' parameter vectors.
#'
#' @param expr     An expression, for example a function body.
#' @param out      Boolean vector storing the output indicators. Initialized to \code{FALSE} and this
#'   should not be changed.
#' 
#' @return Boolean vector of the same format as \code{out}, where \code{TRUE} means the corresponding element
#'   is present in \code{expr}, and \code{FALSE} means it is not.
#'   The elements are \code{y.name}, \code{y.nr} and \code{y.comp} indicating whether the vector \code{y}
#'   is indexed by names (strings), numbers or by more complicated expressions (e.g. arithmetic expressions
#'   or functions). The element \code{y.vec} indicates whether \code{y} is referred to as a vector.
#'   The element \code{y.brc} indicates whether elements of \code{y} are referred to with a single bracket ('[]') or '$'.
#'   More than one of these may be true, if multiple reference methods are used in \code{expr}.
#'   The elements \code{p.name}, \code{p.nr}, \code{p.comp} and \code{p.vec} do the same for \code{p}, and likewise
#'   for \code{theta}, \code{eta} and \code{eps}.
#'   Elements \code{with}, \code{return} and \code{FU} indicate the presence of the associated functions.
#'   If \code{expr} contains illegal components (anything other than a call, pairlist, atomic or name),
#'   then execution is stopped.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
analyzeBody <- function(expr, out = c(y.name = FALSE, y.nr = FALSE, y.comp = FALSE, y.vec = FALSE, y.brc = FALSE,
                                      p.name = FALSE, p.nr = FALSE, p.comp = FALSE, p.vec = FALSE, p.brc = FALSE,
                                      theta.name = FALSE, theta.nr = FALSE, theta.comp = FALSE, theta.vec = FALSE, theta.brc = FALSE,
                                      eta.name = FALSE, eta.nr = FALSE, eta.comp = FALSE, eta.vec = FALSE, eta.brc = FALSE,
                                      eps.name = FALSE, eps.nr = FALSE, eps.comp = FALSE, eps.vec = FALSE, eps.brc = FALSE,
                                      with = FALSE, return = FALSE, FU = FALSE)) {
  varNames <- c("y", "p", "theta", "eta", "eps")  # Variables to check indexing for.
  if (is.name(expr) && as.character(expr) %in% varNames) {
    # y or p appears as vector, i.e., without indexing:
    out[[paste0(as.character(expr), ".vec")]] <- TRUE
    return(out)
  } 
  if (is.atomic(expr) || is.name(expr)) {
    # atomic or name other than y or p, so reached leaf of tree and nothing to change:
    return(out)
  } else if (is.call(expr)) {
    if (length(expr) <= 1) {
      # Call without arguments (probably {}), so nothing to do.
      return(out)
    }
    # Check if expr equals y[[...]] or p[[...]]:
    expr2 <- as.character(expr[[2]])
    if(length(expr) >= 3 && (identical(expr[[1]], quote(`[`)) || identical(expr[[1]], quote(`$`))) && expr2 %in% varNames) {
      out[paste0(expr2, ".brc")] <- TRUE
      return(out)
    }
    if(length(expr) >= 3 && identical(expr[[1]], quote(`[[`))) {
      if(expr2 %in% varNames) {
        cls <- class(expr[[3]])         # index is character string or number or more complicated
        type <- if (cls == "character") "name" else if (cls == "numeric") "nr" else "comp"
        out[[paste0(expr2, ".", type)]] <- TRUE 
        return(out)
      }
    }
    # Check if expr equals "with()" or "return()" or "FU_...()"
    fun <- as.character(expr[[1]])
    if(identical(fun, "with") || identical(fun, "return") || grepl("^FU_", fun[[1]])) {
      if (grepl("^FU_", fun[[1]])) fun <- "FU"
      out[[fun]] <- TRUE
      return(out)
    }
    
    # Otherwise, dig deeper
    outList <- lapply(X = expr, FUN = analyzeBody, out = out)  # Call analyzeBody recursively for other expressions
    return(Reduce("|", x = outList, init = out))
  } else if (is.pairlist(expr)) {
    outList <- unlist(lapply(X = expr, FUN = analyzeBody, pList = out)) # Call analyzeBody recursively for pairlists
    return(Reduce("|", x = outList, init = out))
  } else {
    # User supplied incorrect input
    stop("Don't know how to handle type ", typeof(expr), call. = FALSE)
  }
}


#------------------------------------------ checkVectorUse ------------------------------------------
#' check vector use in a function body
#'
#' Dissects a given function body (technically, any expression) to determine if forbidden vectorizing functions
#' \code{c} or \code{`:`} occur anywhere outside the return value.
#' It does so by going recursively through the expression tree.
#'
#' @param expr            An expression, for example a function body.
#' @param atTailHeadBlock Boolean vector storing whether this expression is at the tail of the function, i.e.,
#'   whether it contains the return statement.
#'   Initialized to \code{TRUE} and this should not be changed.
#' 
#' @return \code{TRUE} if \code{c} or \code{`:`} occur anywhere outside the return statement, and \code{FALSE}
#'   if they only occur in the return value (or not at all).
#'   If \code{expr} contains illegal components (anything other than a call, pairlist, atomic or name),
#'   then the return value is \code{NA}. This should never happen.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
checkVectorUse <- function(expr, atTailHeadBlock = TRUE) {
  # atHeadFinalBlock indicates whether this expression is at the tail of the function, i.e., whether it contains the return statement. Should initialize to TRUE for function body.
  if (is.null(expr)) return(FALSE)        # Empty expression. This is ok
  if (is.atomic(expr)) return(FALSE)      # Atomic expression. This is ok
  if (is.name(expr)) return(FALSE)        # Name (i.e., symbol). This is ok
  if(is.call(expr) && identical(expr[[1]], quote(`{`))) {
    # Expression is block start:
    if(length(expr) == 1) return(FALSE)   # Empty block. This is ok
    out <- FALSE
    if(length(expr) > 2) {
      # More than one statement. Check the ones before the last one:
      for (i in 2:(length(expr)-1)) {
        out <- out || checkVectorUse(expr[[i]], FALSE)   # out to TRUE if any of these is a problem
      }
    }
    return(out || checkVectorUse(expr[[length(expr)]], atTailHeadBlock))  # Evaluate final statement of block
  }
  if(is.call(expr) && (identical(expr[[1]], quote(`c`)) || identical(expr[[1]], quote(`:`)))) {
    # Expression is vector
    if(atTailHeadBlock) return(FALSE) else if(length(expr) > 2) return(TRUE) else return(FALSE)
  }
  if(is.call(expr)) {
    # Any other expression. Check its parts:
    out <- FALSE
    if(length(expr) > 1) {
      for (i in 2:(length(expr))) {
        out <- out || checkVectorUse(expr[[i]], FALSE)   # out to TRUE if any of these is a problem
      }
    }
    return(out)
  }
  if(is.pairlist(expr)) return(TRUE)  # Pairlists are not allowed anywhere in the function body.
  return(NA)  # This should not happen
}


#------------------------------------------ checkModelDef ------------------------------------------
#' Check formatting of model components for symbolic manipulation
#'
#' Checks whether model components (\code{model}, \code{p}, \code{init}, \code{output}) 
#' are formatted correctly for symbolic manipulation.
#' The following checks are performed:
#' * Each component has the correct formals.
#' * The components do not contain "with" or "return" statements, or functions starting with "FU_".
#'   The latter are used to represent the state variables in internal calculations and therefore may not be used for other purposes.
#' * The components do not refer to y, p, theta, eta or eps as vector or indexed via "[.]".
#' * The components refer to elements of y, p, theta, eta or eps in the same way, i.e., either by name or number.
#' * The vectors \code{theta}, \code{eta} and \code{eps} are consistent with this (i.e. whether they name their components or not).
#' * The outputs of the functions \code{init} and \code{p} are consistent with this (i.e. whether they name their outputs or not).
#'
#' @param model    Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
#'   This function should return the numeric vector dy/dt.
#' @param p        Function(theta, eta) of population and individual parameters specifying the model parameters.
#'   This function should return a named numeric vector.
#' @param init     Function(p) of parameters specifying the initial state (at time 0) of the model.
#' @param output   Function(y, p, eps) of state variables, parameters and residual errors specifying the model outputs.
#'   This function should return a numeric vector.
#' @param theta    Population parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param eta      Individual parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param eps      Residual parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param showWarn \code{TRUE} if warnings about missing random parameters are to be shown, \code{FALSE} if not.
#'                 The latter is of use in case only the structural model is of interest.
#'
#' @note \code{theta}, \code{eta} and \code{eps} are only used to generate some output values of the functions
#'   \code{p} and \code{init}, so that their naming can be checked, and to check whether they have names or not.
#'   Their specific values usually do not matter (but see next note).
#' @note These checks are meant to catch the most common mistakes but are not completely failsafe. For example, 
#'   to check consistency in naming they determine if \code{p} and \code{init} provide named or unnamed return values for
#'   the given \code{theta} and \code{eta}. This may fail to detect errors in case the naming depends on the parameter values.
#'   Another example is that the validity or consistency of names or indices is not checked
#'   (e.g., one model component may use different names than another, or indices can be negative).
#' 
#' @return \code{TRUE} if the model components pass all tests, \code{FALSE} if not.
#'   In the latter case, error messages will be displayed.
#'   All encountered errors are shown, but if the formals are incorrect, the other checks are not performed.
#'   In the first case, warning messages may be displayed.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
checkModelDef <- function(model, p, init, output, theta, eta, eps, showWarn) {
  varNames <- c("y", "p", "theta", "eta", "eps")  # Variables to check indexing for.
  errors    <- c()
  warns     <- c()
  fus  <- list(model = model, init = init, p = p, output = output)
  # Check formals:
  frms <- lapply(fus, function(fu) names(formals(fu)))
  expFrms <- list(model = c("t", "y", "p"), init = c("p"), p = c("theta", "eta"), output = c("y", "p", "eps"))
  for (nam in names(frms)) {
    if (!identical(frms[[nam]], expFrms[[nam]])) {
      errors <- c(errors, paste0("\t'", nam, "' function arguments are (", paste0(frms[[nam]], collapse = ", "),
                                 "), but should be (", paste0(expFrms[[nam]], collapse = ", "), ").\n"))
    }
  }
  if (length(errors) == 0) {
    # Check body:
    bodychecks <- lapply(fus, function(fu) analyzeBody(body(fu)))
    bodychecks <- Reduce("|", bodychecks)  # Combine checks of the functions
    # Check special functions:
    for (elt in c("with", "return", "FU")) {
      if (bodychecks[[elt]]) {
        # special function appears:
        eltNam <- if(elt == "FU") "FU_..." else elt
        errors <- c(errors, paste0("\tfunction '", eltNam, "' appears in function body. This is not allowed.\n"))
      }
    }
    # Check indexing of variables (y, p, theta, eta, eps):
    for (elt in varNames) {
      if (bodychecks[[paste0(elt, ".vec")]]) {
        # y or p appears as vector
        errors <- c(errors, paste0("\targument '", elt, "' appears as vector.\n"))
      }
      if (bodychecks[[paste0(elt, ".brc")]]) {
        # y or p appears as vector
        errors <- c(errors, paste0("\targument '", elt, "' indexed by '[.]' or '$'. May only use '[[.]]'.\n"))
      }
      if (bodychecks[[paste0(elt, ".comp")]]) {
        # composite indexing of y or p
        errors <- c(errors, paste0("\targument '", elt, "' has composite indexing. Index by name or number only.\n"))
      }
      if (bodychecks[[paste0(elt, ".name")]] && bodychecks[[paste0(elt, ".nr")]]) {
        # both name and number indexing of y or p
        errors <- c(errors, paste0("\targument '", elt, "' has indexing by name and number. Index by either one, but not both.\n"))
      }
      if (!bodychecks[[paste0(elt, ".name")]] && !bodychecks[[paste0(elt, ".nr")]]) {
        # neither name nor number indexing of y or p. This is not an error
        warns <- c(warns, paste0("\targument '", elt, "' does not appear. Not an error, but you may want to double check.\n"))
      }
    }
    
    # Check that function bodies do not contain vectorization (i.e., 'c()' or ':'), except possibly in the return value:
    vectorUse <- lapply(fus, function(fu) checkVectorUse(body(fu)))
    for (fu in names(vectorUse)) {
      if(vectorUse[[fu]]) errors <- c(errors, paste0("\tfunction '", fu, "' uses vectorization ('c()' or ':'). This may only be used in the return value.\n"))
    }
    
    # Check names of output of "init" and "p" functions:
    selFus <- c("y" = "init", "p" = "p")
    hasOutNames <- {
      # Determine whether the output of functions init and p uses names (TRUE) or not (FALSE):
      parms <- try(p(theta, eta), silent = TRUE)
      vars <- try(init(parms), silent = TRUE)
      # "try" is necessary in case init uses names and p does not provide them. Note that this mistake is caught below for elt = 'p' in case of vars, or elt = 'theta' or 'eta' in case of parms.
      c("y" = if(inherits(vars, "try-error")) NA else !is.null(names(vars)),
        "p" = if(inherits(parms, "try-error")) NA else !is.null(names(parms)))
    }
    #    # Alternative, using function body, assuming final line is the return value:
    #    hasOutNames <- lapply(fus[selFus], function(fu) {
    #      out <- body(fu)
    #      !is.null(names(out[[length(out)]]))
    #    })
    #    names(hasOutNames) <- names(selFus)
    for(elt in names(selFus)) {
      selFu <- selFus[[elt]]
      if(!is.na(hasOutNames[[elt]])) {
        if(hasOutNames[[elt]] && bodychecks[[paste0(elt, ".nr")]]) {
          # names in function output but numbers used in body
          errors <- c(errors, paste0("\targument '", elt, "' indexes by number, but has named output in the '", selFu,
                                     "' function. This needs to be consistent.\n"))
        }
        if(!hasOutNames[[elt]] && bodychecks[[paste0(elt, ".name")]]) {
          # names not in function output but used in body
          errors <- c(errors, paste0("\targument '", elt, "' indexes by name, but has no named output in the '", selFu,
                                     "' function. This needs to be consistent.\n"))
        }
      }
    }
    # Check names of theta, eta, eps
    pars <- list("theta" = theta, "eta" = eta, "eps" = eps)
    hasOutNames <- unlist(lapply(pars, function(par) {
      if (is.null(par)) NA else !is.null(names(par))
    }))
    for (elt in names(hasOutNames)) {
      if(!is.na(hasOutNames[[elt]])) {
        if (hasOutNames[[elt]] && bodychecks[[paste0(elt, ".nr")]]) {
          # names in parameter vector but numbers used in body
          errors <- c(errors, paste0("\tparameter vector '", elt, "' indexes by number, but has named components. This needs to be consistent.\n"))
        }
        if(!hasOutNames[[elt]] && bodychecks[[paste0(elt, ".name")]]) {
          # names not in parameter vector but used in body
          errors <- c(errors, paste0("\tparameter vector '", elt, "' indexes by name, but has no named components. This needs to be consistent.\n"))
        }
      }
    }
  }
  if(length(errors) > 0) {
    processErrors(paste0("VariationalEq::checkModelDef:\n", paste0(errors, collapse = ""), "Exiting.\n"))
    return(FALSE)
  }
  if (length(warns) > 0 && showWarn) {
    # print only if no errors and !showWarn:
    processWarns(paste0("VariationalEq::checkModelDef:\n", paste0(warns, collapse = ""), "Continuing.\n"))
  }
  return(TRUE)
}


#------------------------------------------ getNamesInFunction ------------------------------------------
#' Get names that appear on the left hand side in the given function
#'
#' Gets all strings that appear as names in the given function, e.g., as \code{c(name = object)}.
#'
#' @param fu    Function to analyze.
#'
#' @return A character vector of names.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
getNamesInFunction <- function(fu) {
  out <- doGetNames(body(fu))
  sort(setdiff(out, ""))  # Remove empty names and duplicates
}


#------------------------------------------ doGetNames ------------------------------------------
#' workhorse of \link{getNamesInFunction}
#'
#' Gets the strings that appear as names in the given expression.
#'
#' @param expr    An expression, for example the body of a function.
#'
#' @return A character vector of names.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
doGetNames <- function(expr) {
  if (is.atomic(expr) || is.name(expr)) return(names(expr))  # Return names if atomic or name
  c(names(expr), unlist(lapply(expr, doGetNames)))           # Otherwise, search recursively
}


#------------------------------------------ getRefsInFunction ------------------------------------------
#' Get names that appear on the right hand side in the given function
#'
#' Gets all strings that appear as values in the given function, e.g., as \code{x <- "string"}.
#'
#' @param fu    Function to analyze.
#'
#' @return A character vector of values.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
getRefsInFunction <- function(fu) {
  out <- doGetRefs(body(fu))
  sort(setdiff(out, ""))  # Remove empty referred names and duplicates
}


#------------------------------------------ doGetRefs ------------------------------------------
#' workhorse of \link{getRefsInFunction}
#'
#' Gets the strings that appear as names in the given expression.
#'
#' @param expr    An expression, for example the body of a function.
#'
#' @return A character vector of names.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
doGetRefs <- function(expr) {
  # Return name if atomic character.
  if (is.atomic(expr)) return(if (is.character(expr)) as.character(expr) else NULL)
  if (is.name(expr)) return(NULL)
  unlist(lapply(expr, doGetRefs)) # Otherwise, search recursively
}


#------------------------------------------ createNewNamesRefs ------------------------------------------
#' Generate unique new names and references for the given ones not ending in an alphabetic character
#'
#' Finds names not ending in an alphabetic character, and generates unique names that do end in an alphabetic
#' character by repeatedly appending "X".
#'
#' @param oldNames    Character vector of names to be modified.
#'
#' @return A list of two elements. The first is a vector of indices in \code{oldNames} where the name was modified.
#'         The second is a character vector with modified names.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
createNewNamesRefs <- function(oldNames) {
  inds <- grep("[^[:alpha:]]$", oldNames)
  newNames <- oldNames
  if (length(inds) > 0) {
    for (ind in inds) {
      nam <- oldNames[[ind]]
      while (nam %in% c(oldNames, newNames)) nam <- paste0(nam, "X")
      newNames[[ind]] <- nam
    }
  }
  list(indices = inds, newNames = newNames)
}


#------------------------------------------ replaceNamesRefs ------------------------------------------
#' Replace names and references in a given function
#'
#' Replace the names in the function \code{fu} according to \code{replace}.
#' The purpose of this is to replace names not ending in an alphabetic character by ones that do.
#'
#' @param fu          Function to modify.
#' @param replace     Named character vector of the form \code{original name = new name}.
#'
#' @return The function \code{fu} with modified names.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
replaceNamesRefs <- function(fu, replace) {
  body(fu) <- doReplaceNamesRefs(body(fu), replace)
  fu
}


#------------------------------------------ doReplaceNamesRefs ------------------------------------------
#' workhorse of \link{replaceNamesRefs}
#'
#' Replace the names in \code{oldNames} by the corresponding ones in \code{newNames}, in the expression \code{expr}.
#' The purpose of this is to replace names not ending in an alphabetic character by ones that do.
#'
#' @param expr        An expression, for example the body of a function.
#' @param replace     Named character vector of the form \code{original name = new name}.
#'
#' @return The expression \code{expr} with modified names.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
doReplaceNamesRefs <- function(expr, replace) {
  if (is.atomic(expr) || is.name(expr)) {
    if(!is.null(names(expr))) expr <- plyr::rename(expr, replace, warn_missing = FALSE)
    if (is.atomic(expr)&& is.character(expr)) expr <- plyr::revalue(as.character(expr), replace, warn_missing = FALSE)
    return(expr)
  } else {
    if (is.list(expr)) {
      expr <- lapply(expr, doReplaceNamesRefs, replace = replace)  # This branch is needed only for back-replacement from Deriv::Deriv results.
    } else if (is.call(expr)) {
      expr <- as.call(unlist(lapply(expr, doReplaceNamesRefs, replace = replace))) # Otherwise, replace recursively
    } else stop("unknown object")
    if(!is.null(names(expr))) expr <- plyr::rename(expr, replace, warn_missing = FALSE)
    return(expr)
  }
}


#------------------------------------------ calcVariationsSymb ------------------------------------------
#' Calculate variations symbolically
#'
#' Calculates the variational matrix of the given model, that is, the first and optionally second order derivatives of
#' the model outputs with respect to the model parameters, at given times.
#'
#' @param model    Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
#'   This function should return the numeric vector dy/dt.
#' @param p        Function(theta, eta) of population and individual parameters specifying the model parameters.
#'   This function should return a named numeric vector.
#' @param init     Function(p) of parameters specifying the initial state (at time 0) of the model.
#' @param output   Function(y, p, eps) of state variables, parameters and residual errors specifying the model outputs.
#'   This function should return a numeric vector.
#' @param times    Numeric vector of times where variations are to be evaluated.
#'   Should all be >= 0, but do not need to include the dose time 0.
#' @param theta    Population parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param eta      Individual parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param eps      Residual parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param vartheta Vector of names of population parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(theta)}, or \code{NULL} for none.
#' @param vareta   Vector of names of individual parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(eta)}, or \code{NULL} for none.
#' @param vareps   Vector of names of residual parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(eps)}, or \code{NULL} for none.
#' @param secOrd   \code{TRUE} if second order derivatives have to be computed, \code{FALSE} if not.
#' @param chkModel \code{TRUE} if it has to be checked whether model components (\code{model}, \code{p}, \code{init}, \code{output}) 
#'                 are formatted correctly for symbolic manipulation, \code{FALSE} if not.
#' @param showWarn \code{TRUE} if warnings about missing random parameters are to be shown, \code{FALSE} if not.
#'                 The latter is of use in case only the structural model is of interest.
#' @param ...    Named arguments to be passed to \code{\link[deSolve]{lsoda}}.
#'   Can be used for example to pass events or error tolerances.
#'
#' @return A data frame with columns 't' for time, 'i' for the index of the output element, 'y' for the outputs,
#'   'dy_d<v1>' for the first derivatives and (if \code{secOrd == TRUE}) 'd2y_d<v1>_<v2>' for the second derivatives.
#'   In the column names, variables v1 and v2 are replaced by the names in \code{vartheta}, \code{vareta} and \code{vareps}.
#'   The function displays an error and returns \code{NULL} if \code{vartheta} contains elements not in \code{names(theta)},
#'   and likewise for \code{eta} and \code{eps}.
#'   The same happens if the model components are not formatted correctly for symbolic manipulation and \code{chkModel = TRUE},
#'   or if integration of the model fails.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
calcVariationsSymb <- function(model, p, init, output, times, theta, eta, eps,
                               vartheta, vareta, vareps, secOrd, chkModel, showWarn, ...) {
  if (chkModel) {
    if (!checkModelDef(model, p, init, output, theta, eta, eps, showWarn)) return(NULL)
  }
  # Adapt names appearing in the functions, so that all end in an alphabetic character (this is needed to avoid ambiguity in the results of Deriv::Deriv):
  # First create a list of all names used on the left or right hand side in model components:
  oldNames <- sort(unique(unlist(lapply(list(p, init, model, output), function(fu) {
    c(getNamesInFunction(fu), getRefsInFunction(fu))
  }))))
  newIndsNames <- createNewNamesRefs(oldNames)
  inds <- newIndsNames[[1]]
  newNames <- newIndsNames[[2]]
  # Then, create versions of model components where all names end in an alphabetic character:
  repl <- newNames[inds]
  names(repl) <- oldNames[inds]
  pmdf        <- replaceNamesRefs(p,      repl)
  initmdf     <- replaceNamesRefs(init,   repl)
  modelmdf    <- replaceNamesRefs(model,  repl)
  outputmdf   <- replaceNamesRefs(output, repl)
  thetamdf    <- plyr::rename(theta, repl, warn_missing = FALSE)
  etamdf      <- plyr::rename(eta, repl, warn_missing = FALSE)
  epsmdf      <- plyr::rename(eps, repl, warn_missing = FALSE)
  varthetamdf <- plyr::mapvalues(vartheta, names(repl), repl, warn_missing = FALSE)
  varetamdf   <- plyr::mapvalues(vareta,   names(repl), repl, warn_missing = FALSE)
  varepsmdf   <- plyr::mapvalues(vareps,   names(repl), repl, warn_missing = FALSE)

  # Find parameter values, and times:
  parms <- pmdf(thetamdf, etamdf)
  xTimes <- sort(unique(c(0, times)))
  
  # Create variational function to integrate:
  dfu <- createDerivsSymb(modelmdf, pmdf, initmdf, outputmdf, thetamdf, etamdf, epsmdf,
                          varthetamdf, varetamdf, varepsmdf, secOrd)
  nOutputs <- length(dfu[["output"]](dfu[["init"]](parms), parms, eps))
  fu <- function(t, y, p, eps) {
    c(list(dfu[["model"]](t, y, p)), dfu[["output"]](y, p, eps))
  }
  
  # Integrate variational ODE:
  out <- try(deSolve::lsoda(y = dfu[["init"]](parms), times = xTimes, func = fu, parms = parms, eps = eps, ...))
  if(inherits(out, "try-error")) return(NULL)
  out <- as.data.frame(out[out[, "time"] %in% times, c(1, (ncol(out)-nOutputs+1):ncol(out)), drop = FALSE])
  
  out <- calcVarSymb(out[, -1], dfu[["p"]](thetamdf, etamdf), varthetamdf, varetamdf, varepsmdf, secOrd) # Change to derivs of base parameters (theta, eta, eps)
  colnames <- gsub("\\..*?$", "", names(out))
  colnames <- colnames[!duplicated(colnames)]
  out <- reshape2::melt(cbind(data.frame(t = times), out), id.vars = "t")  # Tall format
  out[, "i"] <- as.numeric(gsub("^.*\\.", "", out[, "variable"]))
  out[, "var"] <- gsub("\\..*?$", "", out[, "variable"])
  out <- reshape2::dcast(out, t+i ~ var, value.var = "value")[, c("t", "i", colnames)]  # Back to wide format
  
  # Revert variable names back to their original values:
  for (ind in inds) {
    colnames(out) <- gsub(paste0("^dy_d\\Q", newNames[[ind]], "\\E$"),
                          paste0("dy_d", oldNames[[ind]]),
                          colnames(out))
    colnames(out) <- gsub(paste0("^d2y_d\\Q", newNames[[ind]], "\\E_(.*)$"),
                          paste0("d2y_d", oldNames[[ind]], "_\\1"),
                          colnames(out))
    colnames(out) <- gsub(paste0("^d2y_d(.*)_\\Q", newNames[[ind]], "\\E$"),
                          paste0("d2y_d\\1_", oldNames[[ind]]),
                          colnames(out))
  }
  
  return(out)
}


#------------------------------------------ isValidVariations ------------------------------------------
#' Check validity of a variational matrix
#'
#' Checks whether the provided input is a valid variational matrix,
#' that is, whether it is in the format as produced by \code{\link{calcVariations}}, \code{\link{calcFirstVariations}}
#' or \code{\link{calcVariationsFim}}.
#' In particular, it is checked whether the input is a data frame of numeric values, with attributes 'type', 'symbolic',
#' 'secOrd', where 'type' is set to 'VariationalMatrix' and 'symbolic' and 'secOrd' are booleans, and with columns
#' 't', 'i', 'y', 'dy_d<v1>' and (if second order) 'd2y_d<v1>_<v2>', where variables v1 and v2 are replaced by names
#' that need to be present in 'theta', 'eta' or 'eps'.
#' The consistency of attributes 'theta', 'eta', 'eps' and 'parameterNormalized' is checked.
#' It is not checked whether attributes 'theta', 'eta', 'eps', 'parameterNormalized' and 'outputNormalized' are present,
#' because these parameters may be NULL in which case they are not listed.
#'
#' @param df     Any object.
#' @param nmdf   Name (string) by which this object should be referred in any error messages. By default "input 'df'".
#'
#' @return \code{NULL} if \code{df} is a valid data frame with a variational matrix,
#'   and an error message or list of error messages if not.
#'
#' @author Martijn van Noort
#' 
#' @noRd 
isValidVariations <- function(df, nmdf = "input 'df'") {
  out <- NULL
  if(!"keepattr" %in% class(df)) {
    return(paste0(nmdf, " is not of the 'keepattr' class."))
  }
  if(!all(c(typeAttr, symbAttr, secOrdAttr) %in% names(attributes(df))) || attr(df, typeAttr) != varTypeVar ||
     !is.logical(attr(df, symbAttr)) || !is.logical(attr(df, secOrdAttr))) {
    return(paste0(nmdf, " does not have the correct attributes of a variational matrix."))
  }
  if (!is.data.frame(df)) {
    return(paste0(nmdf, " should be a data frame."))
  }
  if (!all(unlist(lapply(df, is.numeric)))) {
    out <- c(out, list(paste0(nmdf, " contains non-numeric values.")))
  }
  for (var in c("t", "i", "y")) {
    if (!var %in% names(df)) {
      out <- c(out, list(paste0(nmdf, " does not contain column '", var, "'.")))
    }
  }
  varcols <- grep("^dy_", names(df), value = TRUE)
  if (length(varcols) < 1) {
    out <- c(out, list(paste0(nmdf, " does not contain any columns dy_d<v1>.")))
  }
  secOrd <- attr(df, secOrdAttr)
  if (secOrd) varcols <- c(varcols, grep("^d2y_", names(df), value = TRUE))
  if (length(varcols) != length(unique(varcols))) {
    out <- c(out, list(paste0(nmdf, " contains duplicate columns dy_d<v1> or d2y_d<v1>_<v2>.")))
  }
  if (length(varcols) != ncol(df) - 3) {
    out <- c(out, list(paste0(nmdf, " contains other columns than t, i, y",
                              if (secOrd) ", dy_d<v1> or d2y_d<v1>_<v2>." else " or dy_d<v1>.")))
  }
  nams <- c(names(attr(df, thetaAttr)), names(attr(df, etaAttr)), names(attr(df, epsAttr)))
  namsNorm <- names(attr(df, parNormAttr))
  if((!is.null(nams) || !is.null(namsNorm)) && (length(nams) != length(namsNorm) || !all(nams %in% namsNorm) || !all(namsNorm %in% nams))) {
    out <- c(out, list(paste0(nmdf, " has mismatch between parameter and normalization attributes.")))
  }
  expectvarcols <- if(is.null(nams)) NULL else paste0("dy_d", nams)
  if (secOrd) expectvarcols <- c(expectvarcols, if(is.null(nams)) NULL else paste0("d2y_d", nams, "_", rep(nams, each = length(nams))))
  if (!all(varcols %in% expectvarcols)) {
    out <- c(out, list(paste0(nmdf, " contains unexpected columns ",
                              paste0(setdiff(varcols, expectvarcols), collapse = ", "), ".")))
  }
  return(out)
}


#--------------------------------------------------------------------------------------------------------
#------------------------------------------ Exported functions ------------------------------------------
#--------------------------------------------------------------------------------------------------------


#------------------------------------------ calcVariations ------------------------------------------
#' Calculate variations
#'
#' Calculates the variational matrix of the given model, that is, the first and optionally second order derivatives of
#' the model outputs with respect to the model parameters, at given times.
#' The user can choose whether the computation uses numerical or symbolic derivatives.
#'
#' @param model       Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
#'   This function should return the numeric vector dy/dt.
#'    See the notes for requirements in case of symbolic derivations.
#' @param p           Function(theta, eta) of population and individual parameters specifying the model parameters.
#'   This function should return a named numeric vector.
#'    See the notes for requirements in case of symbolic derivations.
#' @param init        Function(p) of parameters specifying the initial state (at time 0) of the model.
#'    See the notes for requirements in case of symbolic derivations.
#' @param output      Function(y, p, eps) of state variables, parameters and residual errors specifying the model outputs.
#'   This function should return a numeric vector.
#'    See the notes for requirements in case of symbolic derivations.
#' @param times       Numeric vector of times where variations are to be evaluated.
#'   Should all be >= 0, but do not need to include the dose time 0.
#' @param theta       Population parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param eta         Individual parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param eps         Residual parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param vartheta    Vector of names of population parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(theta)}, or \code{NULL} for none.
#'   By default equal to \code{names(theta)}.
#' @param vareta      Vector of names of individual parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(eta)}, or \code{NULL} for none.
#'   By default equal to \code{names(eta)}.
#' @param vareps      Vector of names of residual parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(eps)}, or \code{NULL} for none.
#'   By default equal to \code{names(eps)}.
#' @param secOrd      \code{TRUE} if second order derivatives have to be computed, \code{FALSE} (default) if not.
#' @param symbolic    \code{TRUE} (default) if derivatives are to be computed symbolically, \code{FALSE} if numerically.
#' @param chkModel    \code{TRUE} (default) if it has to be checked whether model components (\code{model}, \code{p}, \code{init}, \code{output}) 
#'   are formatted correctly for symbolic manipulation, \code{FALSE} if not.
#' @param showWarn    \code{TRUE} (default) if warnings about missing random parameters are to be shown, \code{FALSE} if not.
#'   The latter is of use in case only the structural model is of interest.
#' @param ...         Named arguments to be passed to \code{\link[deSolve]{lsoda}}.
#'   Can be used for example to pass events or error tolerances.
#'
#' @return A data frame with columns 't' for time, 'i' for the index of the output element, 'y' for the outputs,
#'   'dy/d<v1>' for the first derivatives and (if \code{secOrd == TRUE}) 'd2y/d<v1>_<v2>' for the second derivatives.
#'   In the column names, variables v1 and v2 are replaced by the names in \code{vartheta}, \code{vareta} and \code{vareps}.
#'   The data frame has attributes 'theta', 'eta' and 'eps' listing the variables \code{theta}, \code{eta} and \code{eps},
#'   as named vectors.
#'   The 'type' attribute is set to 'VariationalMatrix', and the 'secOrd' and 'symbolic' ones to the values of
#'   \code{secOrd} and \code{symbolic}, respectively.
#'   The variational matrix is not normalized, as reflected in the attributes 'parameterNormalized' and 'outputNormalized',
#'   that record input and output normalization, respectively.
#'   The attribute 'parameterNormalized' is set to a named boolean vector of \code{FALSE}, with the theta, eta and eps
#'   parameter names as names.
#'   The attribute 'outputNormalized' is set to \code{FALSE}.
#'   The function displays an error and returns \code{NULL} if \code{vartheta} contains elements not in \code{names(theta)},
#'   and likewise for \code{eta} and \code{eps}, or if \code{times} includes a negative value.
#'
#' @details If \code{symbolic==TRUE}, then the variational matrix is calculated by first taking symbolic derivatives of
#'   the differential equations and associated functions (\code{p}, \code{init} and \code{output}), and solving the resulting
#'   system of differential equations numerically.
#'   If \code{symbolic==FALSE}, then the different equations are solved numerically, and subsequently their derivatives are
#'   computed numerically (using Richardson's method, \url{https://en.wikipedia.org/wiki/Richardson_extrapolation}).
#'   The first method tends to give more stable results, and hence is set as default.
#'   It is not clear how the complexities of these methods compare. The first requires solving a system of differential
#'   equations of higher dimension, while the second requires multiple solutions to the original system (of lower dimension).
#' 
#' @note It is recommended to use symbolic derivation (\code{symbolic==TRUE}) rather than numerical derivation,
#'   as it is more accurate, especially for derivatives that are formally zero.
#'   In numerical calculation, their value may be estimated as close to zero rather than exactly zero,
#'   and this may give the wrong answers, especially for the sensitivity and FIM methods.
#'   Downsides of symbolic derivation are that it can be considerably slower than the numerical option,
#'   and it places more constraints on the definition of the functions \code{model}, \code{p},
#'   \code{init} and \code{output}.
#'   Therefore it is advisable to compare symbolic and numerical versions,
#'   in order to spot any mistakes in these definitions.
#'   The constraints are necessary so that derivatives are properly computed.
#'   They are checked (to an extent) in case \code{chkModel = TRUE}.
#'   They are:
#'   * The function arguments of these functions should be named exactly as specified in the arguments above.
#'     I.e., it is wrong to define \code{model <- function(t, x, q) ...}.
#'   * Do not use \code{with} in the function definitions.
#'   * Do not use a \code{return} statement in the function definitions.
#'     It is allowed to use composite function bodies (i.e., containing multiple statements).
#'     The final statement is the return value.
#'   * Refer to function arguments only component-wise, i.e., as \code{y[["a"]]}, \code{y[[1]]} etc, and not as vectors.
#'     A function definition of the form \code{p <- function(theta, eta) theta} should not be used.
#'     Rather use \code{p <- function(theta, eta) c(theta[["CL"]], theta[["V"]])}.
#'     Components may be referred to by \code{[[]]}, but not by \code{[]} or \code{$}.
#'   * Components of the state vector should be referred to either by index or by name, but not both, also not between functions.
#'     Likewise for components of the parameter vectors.
#'   * Do not define any functions with names starting with "FU_". These are used internally to represent state
#'     variables and therefore may not be used for other purposes.
#' 
#' @note The functions \code{model}, \code{p}, \code{init} and \code{output} should satisfy the rules imposed by
#'   \code{\link[deSolve]{lsoda}}, where \code{output} corresponds to the additional output of the \code{lsoda} function \code{func}.
#' 
#' @note If this function is used for sensitivity or aliasing calculations, then individual and random parameters are not used,
#'   so the function \code{p} and \code{output} may ignore these inputs, and \code{eta} and \code{eps} may be set to \code{NULL}.
#'   Furthermore, in this case \code{secOrd = FALSE}. See \code{\link{calcFirstVariations}}.
#'
#' @note If this function is used to calculate the Fisher Information matrix, then derivatives are evaluated at \code{eta = 0}
#'   and \code{eps = 0}, so the values of these inputs should be set to 0. Furthermore, in this case \code{secOrd = TRUE}.
#'   See \code{\link{calcVariationsFim}}.
#'
#' @export
#' 
#' @family calculation functions
#'
#' @author Martijn van Noort
calcVariations <- function(model, p, init, output, times, theta, eta, eps,
                           vartheta = names(theta), vareta = names(eta), vareps = names(eps), secOrd = FALSE,
                           symbolic = TRUE, chkModel = TRUE, showWarn = TRUE, ...) {
  times <- unique(times)
  # Check that vartheta, vareta and vareps contain valid parameter names:
  for (nam in c("theta", "eta", "eps")) {
    if (length(setdiff(get(paste0("var", nam)), names(get(nam)))) > 0) {
      processErrors(paste0("VariationalEq::calcVariations: 'var", nam, "' should be a subset of (or equal to) '", nam, "'.\nExiting.\n"))
      return(NULL)
    }
  }
  if (any(times < 0)) {
    tms <- sort(times[times < 0])
    dts <- if (length(tms) > 4) ", ..." else ""
    processErrors(paste0("VariationalEq::calcVariations: 'times' contains negative (i.e., pre-initialization) times ",
                         paste0(head(tms, 4), collapse = ", "), dts, ".\nExiting.\n"))
    return(NULL)
  }
  out <- if (symbolic) {
    calcVariationsSymb(model, p, init, output, times, theta, eta, eps, vartheta, vareta, vareps, secOrd, chkModel, showWarn, ...)
  } else {
    calcVariationsNum(model, p, init, output, times, theta, eta, eps, vartheta, vareta, vareps, secOrd, ...)
  }
  if (is.null(out)) return(out)
  out <- out[order(out[, "i"], out[, "t"]), ]
  attr(out, typeAttr)    <- varTypeVar
  attr(out, symbAttr)    <- symbolic
  attr(out, secOrdAttr)  <- secOrd
  attr(out, thetaAttr)   <- theta
  attr(out, etaAttr)     <- eta
  attr(out, epsAttr)     <- eps
  nams <- c(names(theta), names(eta), names(eps))
  parNorm <- rep(FALSE, length(nams))
  names(parNorm) <- nams
  attr(out, parNormAttr) <- parNorm
  attr(out, outNormAttr) <- FALSE
  class(out) <- c('keepattr', class(out))   # This keeps attributes when subsetting
  return(out)
}


#------------------------------------------ calcFirstVariations ------------------------------------------
#' Calculate first order variations
#'
#' Calculates the first variational matrix of the given model, that is, the first order derivatives of
#' the model outputs with respect to the model parameters, at given times.
#' This is tailored for use with the sensitivity and aliasing methods.
#' Therefore model parameters are assumed to be given directly, rather than as a combination of population and
#' individual parameters, and residual error parameters are not included.
#'
#' @param model       Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
#'   This function should return the numeric vector dy/dt.
#' @param parms       Parameter values, as named numeric vector. May be \code{NULL} if none are needed
#'   (but then no derivatives will be calculated).
#' @param init        Function(p) of parameters specifying the initial state (at time 0) of the model.
#' @param outputPred  Function(y, p) of state variables and parameters specifying the model prediction, i.e., the
#'   output without residual.
#'   This function should return a numeric vector.
#' @param times       Numeric vector of times where variations are to be evaluated.
#'   Should all be >= 0, but do not need to include the dose time 0.
#' @param varp        Vector of names of parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(parms)}, or \code{NULL} for none.
#'   By default equal to \code{names(parms)}.
#' @param symbolic    \code{TRUE} (default) if derivatives are to be computed symbolically, \code{FALSE} if numerically.
#'   See \code{\link{calcVariations}} for details.
#' @param chkModel    \code{TRUE} (default) if it has to be checked whether model components (\code{model}, \code{parms}, \code{init}, \code{outputPred}) 
#'   are formatted correctly for symbolic manipulation, \code{FALSE} if not.
#' @param ...         Named arguments to be passed to \code{\link[deSolve]{lsoda}}.
#'   Can be used for example to pass events or error tolerances.
#'
#' @return A data frame with columns 't' for time, 'i' for the index of the output element, 'y' for the outputs,
#'   and 'dy/d<v1>' for the first derivatives.
#'   In the column names, variable v1 is replaced by the names in \code{varp}.
#'   The data frame has an attribute 'theta' listing the variables in \code{parms} as a named vector,
#'   and 'eta' and 'eps' that are empty.
#'   The 'type' attribute is set to 'VariationalMatrix', 'secOrd' to \code{FALSE}, and 'symbolic' to the value of
#'   \code{symbolic}.
#'   The variational matrix is not normalized, as reflected in the attributes 'parameterNormalized' and 'outputNormalized',
#'   that record input and output normalization, respectively.
#'   The attribute 'parameterNormalized' is set to a named boolean vector of \code{FALSE}, with the \code{parms}
#'   names as names.
#'   The attribute 'outputNormalized' is set to \code{FALSE}.
#'   The function displays an error and returns \code{NULL} if \code{varp} contains elements not in \code{names(p)},
#'   or if \code{times} includes a negative value.
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
calcFirstVariations <- function(model, parms, init, outputPred, times, varp = names(parms), symbolic = TRUE, chkModel = TRUE, ...) {
  if (symbolic) {
    # For symbolic calculations, need to define functions for p and output in an explicit way, i.e. without
    # using vector evaluation or reference to other functions:
    thetaComps <- names(parms)
    thetaNames <- TRUE
    if (is.null(thetaComps)) {
      thetaComps <- seq_along(theta) 
      thetaNames <- FALSE
    }
    # thetaComps contains names or indices of theta components
    argVec <- if(thetaNames) paste0("\"", thetaComps, "\" = theta[[\"", thetaComps, "\"]]") else paste0("theta[[", thetaComps, "]]")
    argVec <- paste0("c(", paste0(argVec, collapse = ", "), ")")
    argVec <- parse(text = argVec)[[1]]
    pfunc <- function(theta, eta) NULL
    body(pfunc) <- argVec
    outfunc <- function() NULL
    formals(outfunc) <- c(formals(outputPred), list(eps = expr()))
    body(outfunc) <- body(outputPred)
    calcVariations(model = model, p = pfunc, init = init, output = outfunc,
                   times = times, theta = parms, eta = NULL, eps = NULL, vartheta = varp, vareta = NULL,
                   vareps = NULL, secOrd = FALSE, symbolic = symbolic,
                   chkModel = chkModel, showWarn = FALSE, ...)
  } else {
    # For numerical calculations, can use a simpler approach:
    calcVariations(model = model, p = function(theta, eta) { theta }, init = init,
                   output = function(y, p, eps) { outputPred(y, p) },
                   times = times, theta = parms, eta = NULL, eps = NULL, vartheta = varp, vareta = NULL,
                   vareps = NULL, secOrd = FALSE, symbolic = symbolic,
                   chkModel = chkModel, showWarn = FALSE, ...)
  }
}


#------------------------------------------ calcVariationsFim ------------------------------------------
#' Calculate first and second order variations for FIM
#'
#' Calculates the variational matrix of the given model, that is, the first and second order derivatives of
#' the model outputs with respect to the model parameters, at given times.
#' This is tailored for use with the FIM methods.
#' Therefore derivatives of individual and residual error parameters are evaluated at 0, and so only their
#' names need to be provided as input.
#'
#' @param model       Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
#'   This function should return the numeric vector dy/dt.
#' @param p           Function(theta, eta) of population and individual parameters specifying the model parameters.
#'   This function should return a named numeric vector.
#' @param init        Function(p) of parameters specifying the initial state (at time 0) of the model.
#' @param output      Function(y, p, eps) of state variables, parameters and residual errors specifying the model outputs.
#'   This function should return a numeric vector.
#' @param times       Numeric vector of times where variations are to be evaluated.
#'   Should all be >= 0, but do not need to include the dose time 0.
#' @param theta       Population parameter values, as named numeric vector. May be \code{NULL} if none are needed.
#' @param nmeta       Names of individual parameters, as string vector. May be \code{NULL} if none are needed.
#' @param nmeps       Names of residual parameter values, as string vector. May be \code{NULL} if none are needed.
#' @param vartheta    Vector of names of population parameters for which variations are to be calculated.
#'   Should be a subset of \code{names(theta)}, or \code{NULL} for none.
#'   By default equal to \code{names(theta)}.
#' @param vareta      Vector of names of individual parameters for which variations are to be calculated.
#'   Should be a subset of \code{nmeta}, or \code{NULL} for none.
#'   By default equal to \code{nmeta}.
#' @param vareps      Vector of names of residual parameters for which variations are to be calculated.
#'   Should be a subset of \code{nmeps}, or \code{NULL} for none.
#'   By default equal to \code{nmeps}.
#' @param symbolic    \code{TRUE} (default) if derivatives are to be computed symbolically, \code{FALSE} if numerically.
#'   See \code{\link{calcVariations}} for details.
#' @param chkModel    \code{TRUE} (default) if it has to be checked whether model components (\code{model}, \code{p}, \code{init}, \code{output}) 
#'   are formatted correctly for symbolic manipulation, \code{FALSE} if not.
#' @param ...         Named arguments to be passed to \code{\link[deSolve]{lsoda}}.
#'   Can be used for example to pass events or error tolerances.
#'
#' @return A data frame with columns 't' for time, 'i' for the index of the output element, 'y' for the outputs,
#'   'dy/d<v1>' for the first derivatives and 'd2y/d<v1>_<v2>' for the second derivatives.
#'   In the column names, variables v1 and v2 are replaced by the names in \code{vartheta}, \code{vareta} and \code{vareps}.
#'   The data frame has an attribute 'theta' listing the values of \code{theta} as a named vector, and attributes
#'   'eta' and 'eps' containing named vectors with 0 values for all parameters in \code{nmeta} and \code{nmeps}.
#'   The 'type' attribute is set to 'VariationalMatrix', 'secOrd' to \code{TRUE} and 'symbolic' to the value of
#'   \code{symbolic}.
#'   The variational matrix is not normalized, as reflected in the attributes 'parameterNormalized' and 'outputNormalized',
#'   that record input and output normalization, respectively.
#'   The attribute 'parameterNormalized' is set to a named boolean vector of \code{FALSE}, with the theta, eta and eps
#'   parameter names as names.
#'   The attribute 'outputNormalized' is set to \code{FALSE}.
#'   The function displays an error and returns \code{NULL} if \code{vartheta} contains elements not in \code{names(theta)},
#'   and likewise for \code{nmeta} and \code{nmeps}, or if \code{times} includes a negative value.
#'
#' @note The inputs \code{vareta} and \code{vareps} should contain all random parameters with nonzero variance, even if the
#'   FIM calculation considers the variance and covariances of such a parameter as fixed (that is, derivatives with respect to
#'   these variance and covariances do not appear in the FIM).
#'   That is, random parameters may be dropped from \code{vareta} and \code{vareps} only if they have zero variance.
#'   Notice that any variance and covariance parameters appearing in the FIM necessarily correspond to random parameters with
#'   nonzero variance, so inclusion of all random parameters with nonzero variance is not only necessary but also sufficient.
#'   That is, \code{vareta} and \code{vareps} can be set equal to the sets of all individual and residual parameters with
#'   nonzero variance, respectively.
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
calcVariationsFim <- function(model, p, init, output, times, theta, nmeta, nmeps,
                              vartheta = names(theta), vareta = nmeta, vareps = nmeps, symbolic = TRUE,
                              chkModel = TRUE,...) {
  eta <- rep(0, length(nmeta))
  names(eta) <- nmeta
  eps <- rep(0, length(nmeps))
  names(eps) <- nmeps
  calcVariations(model = model, p = p, init = init, output = output, times = times,
                 theta = theta, eta = eta, eps = eps, vartheta = vartheta, vareta = vareta, vareps = vareps,
                 secOrd = TRUE, symbolic = symbolic, chkModel = chkModel,
                 showWarn = TRUE, ...)
}


#------------------------------------------ getAllParsVariations ------------------------------------------
#' Get all parameter values of a variational matrix
#'
#' Gets all parameter values that were used to create a given variational matrix.
#'
#' @param df      A data frame containing a variational matrix, optionally with second order derivatives.
#'
#' @return   A named list of three named numeric vectors containing the theta's, eta's and eps's used to
#'   create the variational matrix.
#'   The list elements are named "theta", "eta" and "eps", respectively.
#'   The vector elements have the names of the parameters.
#'   The function displays an error and returns \code{NULL} if \code{df} is not a valid variational matrix.
#'   
#'
#' @export
#' 
#' @family retrieval functions
#'
#' @author Martijn van Noort
getAllParsVariations <- function(df) {
  # Check input validity:
  valid <- isValidVariations(df)
  if (!is.null(valid)) {
    processErrors(paste0("VariationalEq::getAllParsVariations:\n", paste0(paste0("\t", unlist(valid), "\n"), collapse = ""), "Exiting.\n"))
    return(NULL)
  }
  nams <- c("theta", "eta", "eps")
  out <- lapply(nams, function(var) attr(df, var))
  names(out) <- nams
  return(out)
}


#------------------------------------------ getParsVariations ------------------------------------------
#' Get column parameter values of a variational matrix
#'
#' Gets the parameter values of the columns of a variational matrix that were used to create a given variational matrix,
#' as a list of population, individual and residual parameters.
#' This is the subset of \code{getAllParsVariations} corresponding to the columns present in the variational matrix.
#' The output is the same as that of \code{getParVecVariations} except it is formatted as a list of three vectors.
#'
#' @param df      A data frame containing a variational matrix, optionally with second order derivatives.
#'
#' @return   A named list of three named numeric vectors containing the theta's, eta's and eps's used to
#'   create the variational matrix and present as columns of the matrix, in the order of appearance in the matrix.
#'   The list elements are named "theta", "eta" and "eps", respectively.
#'   The vector elements have the names of the parameters.
#'   The function displays an error and returns \code{NULL} if \code{df} is not a valid variational matrix.
#'   
#'
#' @export
#' 
#' @family retrieval functions
#'
#' @author Martijn van Noort
getParsVariations <- function(df) {
  out <- getAllParsVariations(df)
  if(is.null(out)) return(NULL)
  varcols <- grep("^dy_", names(df), value = TRUE)
  varnms <- unique(gsub("^dy_d", "", varcols))  # variables names
  lapply(out, function(vec) vec[varnms[varnms %in% names(vec)]])
}


#------------------------------------------ getParVecVariations ------------------------------------------
#' Get column parameter values of a variational matrix as a single vector
#'
#' Gets the parameter values of the columns of a variational matrix that were used to create a given variational matrix,
#' as a single vector.
#' This is the subset of \code{getAllParsVariations} corresponding to the columns present in the variational matrix.
#' The output is the same as that of \code{getParVariations} except it is formatted as a single vector.
#'
#' @param df      A data frame containing a variational matrix, optionally with second order derivatives.
#'
#' @return   A named numeric vector containing the theta's, eta's and eps's used to
#'   create the variational matrix and present as columns of the matrix, in the order of appearance in the matrix.
#'   The names are the names of the parameters.
#'   The function displays an error and returns \code{NULL} if \code{df} is not a valid variational matrix.
#'   
#'
#' @export
#' 
#' @family retrieval functions
#'
#' @author Martijn van Noort
getParVecVariations <- function(df) {
  out <- getAllParsVariations(df)
  if(is.null(out)) return(NULL)
  varcols <- grep("^dy_", names(df), value = TRUE)
  varnms <- unique(gsub("^dy_d", "", varcols))  # variables names
  names(out) <- NULL
  out <- unlist(out, recursive = FALSE, use.names = TRUE)
  if(!all(varnms %in% names(out))) {
    # Unknown columns in variational matrix. This should not happen:
    processErrors(paste0("VariationalEq::getParVecVariations: no parameter values for columns ",
                         paste0(setdiff(varnms, names(out)), ", "), ".\nExiting.\n"))
    return(NULL)
  }
  out[varnms]
}


#------------------------------------------ normalizeVariations ------------------------------------------
#' Normalize variations
#'
#' Applies parameter- and/or output-normalization to a given first and (optionally) second order variational matrix.
#' For the first order matrix, parameter-normalization is performed by right multiplication with a matrix
#' P = diag(p), where p is a vector of parameter values (if the parameter is to be normalized) or 1 (if not),
#' and the output-normalization involves left multiplication with a matrix Y^(-1) where Y = diag(y) and y is
#' the vector of output values.
#' For the second order matrix (actually tensor), parameter-normalization involves two suitable multiplications
#' with P, and the output-normalization is done in the same way as for the first order matrix.
#' If the matrix was already normalized one way or the other, then the normalization is adapted to the new settings.
#' Parameter-normalization is selected per parameter, while output-normalization is a single choice. 
#'
#' @param df      A data frame containing a variational matrix, optionally with second order derivatives,
#'   as produced by \code{\link{calcVariations}}, \code{\link{calcFirstVariations}} or
#'   \code{\link{calcVariationsFim}}, that is, with columns 't', 'i', 'y', 'dy_d<v1>' and
#'   (if second order) 'd2y_d<v1>_<v2>', where variables v1 and v2 are replaced by names.
#'   Optionally this matrix may be parameter- or output-normalized.
#' @param parNorm Vector (not necessarily named) of booleans of length equal to the number of parameters used to
#'   construct \code{df}, i.e., those listed in its attributes 'theta', 'eta' and 'eps'.
#'   The value \code{TRUE} means the corresponding individual parameter should be normalized in the FIM,
#'   \code{FALSE} means not.
#'   Instead of a vector, may also provide a single boolean, where \code{TRUE} (default) stands for a vector of
#'   all \code{TRUE}, and \code{FALSE} for a vector of all \code{FALSE}.
#'   May also provide a character vector with the names of the parameters to be normalized.
#'   Invalid names are ignored, with a warning.
#' @param outNorm \code{TRUE} (default) if variations should be output-normalized, \code{FALSE} if not.
#'
#' @return   A data frame in same format, where the columns 'dy_d<v1>' contain the values, normalized as specified.
#'   The function displays an error and returns \code{NULL} if the input is not in the specified format,
#'   \code{parNorm} does not contain information on all parameters, or \code{outNorm = TRUE} and the 'y' column
#'   of \code{df} contains zeroes.
#'   The attributes are as those of \code{df}, with 'parameterNormalized' and 'outputNormalized' adapted according
#'   to the new normalization.
#'
#' @export
#' 
#' @family result modifiers
#'
#' @author Martijn van Noort
normalizeVariations <- function(df, parNorm = TRUE, outNorm = TRUE) {
  # Check input validity:
  valid <- isValidVariations(df)
  if (!is.null(valid)) {
    processErrors(paste0("VariationalEq::normalizeVariations:\n", paste0(paste0("\t", unlist(valid), "\n"), collapse = ""), "Exiting.\n"))
    return(NULL)
  }
  parms <- c(attr(df, thetaAttr), attr(df, etaAttr), attr(df, epsAttr))
  varcols <- grep("^dy_", names(df), value = TRUE)
  varnms <- unique(gsub("^dy_d", "", varcols))  # variables names
  if (outNorm && any(df[, "y"] == 0)) {
    processErrors("VariationalEq::normalizeVariations: column 'y' in input 'df' contains zeroes, so cannot normalize output.\nExiting.\n")
    return(NULL)
  }
  # Check parNorm:
  if (is.logical(parNorm)) {
    if (length(parNorm) == 1) parNorm <- rep(parNorm, length(parms))
  } else {
    # parNorm treated as character vector:
    if (!all(parNorm %in% names(parms))) processWarns("VariationalEq::normalizeVariations: 'parNorm' contains unknown parameters. They will be ignored.\nContinuing.\n")
    parNorm <- names(parms) %in% parNorm
  }
  if(length(parNorm) != length(parms)) {
    processErrors(paste0("VariationalEq::normalizeVariations: 'parNorm' has incorrect format.\nExiting.\n"))
    return(NULL)
  }
  names(parNorm) <- names(parms)
  # Existing normalization:
  existParNorm <- attr(df, parNormAttr)
  if (length(varnms) > 0 && !all(varnms %in% names(existParNorm))) {
    processErrors(paste0("VariationalEq::normalizeVariations: No prior normalization information for parameters ",
                  paste0(setdiff(varnms, names(existParNorm))), " .\nExiting.\n"))
    return(NULL)
  }
  
  # Perform normalization:
  nSecOrd <- 0  # nr of second order cases
  secOrd <- attr(df, secOrdAttr)
  # Change parameter normalization:
  for (var in varnms) {
    fac <- if(parNorm[[var]]) parms[[var]] else 1
    if(existParNorm[[var]] && parms[[var]] > 0) fac <- fac/parms[[var]]  # if parsm[[var]] == 0 then prior normalization would have normalized to 0, in which case the value of 'fac' does not matter.
    df[, paste0("dy_d", var)] <- df[, paste0("dy_d", var), drop = FALSE] * fac
    if (secOrd) {
      for (var2 in varnms) {
        varnm <- paste0("d2y_d", var, "_", var2)
        if (varnm %in% names(df)) {
          nSecOrd <- nSecOrd + 1
          fac2 <- if(parNorm[[var2]]) parms[[var2]] else 1
          if(existParNorm[[var2]] && parms[[var2]] > 0) fac2 <- fac2/parms[[var2]]
          df[, varnm] <- df[, varnm, drop = FALSE] * fac * fac2
        }
      }
    }
  }
  if (nSecOrd > 0 & nSecOrd != length(varnms)*(length(varnms)+1)/2) {
    processWarns("VariationalEq::normalizeVariations: unexpected number of second order derivatives.\nThis should not happen, but I will continue for now.\n")
  }
  
  existOutNorm <- attr(df, outNormAttr)
  if ((existOutNorm && !outNorm) || (!existOutNorm && outNorm)) {
    nams <- c(paste0("dy_d", varnms), grep("^d2y_d", names(df), value = TRUE))  # Names of first and second order columns
    df[, nams] <- df[, nams, drop = FALSE] * if(outNorm) 1 / df[, "y", drop = TRUE] else df[, "y", drop = TRUE]
  }
  attr(df, parNormAttr) <- parNorm
  attr(df, outNormAttr) <- outNorm
  return(df)
}


#------------------------------------------ isNormalizedVariations ------------------------------------------
#' Check whether a variational matrix is normalized
#'
#' Checks whether a variational matrix is normalized.
#' The check involves checking the validity of the input and the value of its attributes \code{parameterNormalized}
#' and \code{outputNormalized}.
#'
#' @param df      A data frame containing a variational matrix, optionally with second order derivatives.
#'
#' @return   A named boolean list stating the parameter- or output-normalization of the matrix.
#'   Function displays an error and returns \code{NULL} if the input is not in the specified format.
#'
#' @export
#' 
#' @family checkers
#'
#' @author Martijn van Noort
isNormalizedVariations <- function(df) {
  # Check input validity:
  valid <- isValidVariations(df)
  if (!is.null(valid)) {
    processErrors(paste0("VariationalEq::isNormalizedVariations:\n", paste0(paste0("\t", unlist(valid), "\n"), collapse = ""), "Exiting.\n"))
    return(NULL)
  }
  list("parameterNormalized" = attr(df, parNormAttr), "outputNormalized" = attr(df, outNormAttr))
}


#------------------------------------------ plotVariations ------------------------------------------
#' Plot variations
#'
#' Creates plots of the model and its variations, for a given variational matrix.
#'
#' @param df         A data frame containing a variational matrix with or without second order derivatives,
#'   as produced by \code{\link{calcVariations}}, \code{\link{calcFirstVariations}} or \code{\link{calcVariationsFim}},
#'   that is, with columns 't', 'i', 'y', 'dy_d<v1>' and optionally 'd2y_d<v1>_<v2>', where variables v1 and v2
#'   are replaced by names.
#' @param outNames   A vector of strings containing the names of the outputs y_i in order.
#'   They are used as labels.
#'   If \code{NULL} (default), then the index \code{i} is used as label.
#' @param plotModel  If \code{TRUE} (default), the model is plotted, if \code{FALSE} then not.
#' @param plotFirst  If \code{TRUE} (default), the first variations are plotted, if \code{FALSE} then not.
#' @param plotSecond If \code{TRUE} (default), the second variations are plotted if present, if \code{FALSE} then not.
#' @param plotIIV    If \code{TRUE}, the variations with respect to IIV variables are plotted if present, if \code{FALSE} (default) then not.
#'   This requires the attribute 'eta' to be present.
#' @param plotRes    If \code{TRUE} (default), the variations with respect to residual variables are plotted if present, if \code{FALSE} then not.
#'   This requires the attribute 'eps' to be present.
#' @param plotPoints If \code{FALSE} (default), only lines are plotted, if \code{TRUE} also plots points.
#'
#' @return   A plot.
#'   In case of error, an error message is printed and the return value is \code{NULL}.
#'
#' @export
#' 
#' @family plotting and printing
#'
#' @author Martijn van Noort
plotVariations <- function(df, outNames = NULL, plotModel = TRUE, plotFirst = TRUE, plotSecond = TRUE, plotIIV = FALSE, plotRes = TRUE, plotPoints = FALSE) {
  # Check input validity:
  valid <- isValidVariations(df)
  if (!is.null(valid)) {
    processErrors(paste0("VariationalEq::plotVariations:\n", paste0(paste0("\t", unlist(valid), "\n"), collapse = ""), "Exiting.\n"))
    return(NULL)
  }
  # Create output names:
  if (!is.null(outNames) & length(outNames) != length(unique(df[, "i"]))) {
    processErrors(paste0("VariationalEq::plotVariations: outNames should be of length ", length(unique(df[, "i"])), ".\nExiting.\n"))
    return(NULL)
  }
  df[, "output"] <- if(is.null(outNames)) {
    factor(df[, "i"], levels = sort(unique(df[, "i"])), ordered = TRUE)
  } else {
    factor(outNames[df[, "i"]], levels = outNames, ordered = TRUE)
  }
  # Identify variables to plot:
  varcols <- if(plotModel) "y" else NULL
  varcols <- c(varcols, if(plotFirst) grep("^dy_", names(df), value = TRUE) else NULL)
  varcols <- c(varcols, if(plotSecond) grep("^d2y_", names(df), value = TRUE) else NULL)
  if (!plotIIV && length(attr(df, etaAttr)) > 0) {
    # Remove all columns containing an IIV parameter:
    grepstr <- paste0("_d?(", paste0(names(attr(df, etaAttr)), collapse = "|"), ")")
    varcols <- grep(grepstr, varcols, value = TRUE, invert = TRUE)
  } 
  if (!plotRes && length(attr(df, epsAttr)) > 0) {
    # Remove all columns containing a residual parameter:
    grepstr <- paste0("_d?(", paste0(names(attr(df, epsAttr)), collapse = "|"), ")")
    varcols <- grep(grepstr, varcols, value = TRUE, invert = TRUE)
  } 
  # Set to tall:
  dfTall <- reshape2::melt(data = df[, c("t", "output", varcols)], id.vars = c("t", "output"))
  pl <- ggplot2::ggplot(dfTall, mapping = aes(x = t, y = value, group = interaction(variable, output), color = output)) +
    geom_line()
  if(plotPoints) pl <- pl + geom_point()
  return(pl + facet_wrap(~variable, scales = "free_y") + ggplot2::theme_bw())
}
