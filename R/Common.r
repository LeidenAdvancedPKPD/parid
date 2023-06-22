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


.onAttach <- function(libname, pkgname) {
  packageStartupMessage(pkgname, ", version ", packageVersion(pkgname), ", (C) 2023 LAP&P Consultants BV. See COPYING for license info")
}

# Attributes for variational matrices and FIMs:
typeAttr      <- "type"                # Attribute where the type is stored.
thetaAttr     <- "theta"               # Attribute where theta parameter values of variational matrix or FIM are stored. The attribute contains a named numeric vector with the theta values. The names are the parameter names.
symbAttr      <- "symbolic"            # Attribute storing a boolean whether the calculation was performed symbolically or numerically.

varTypeVar    <- "VariationalMatrix"   # Value of "type" attribute of variational matrix
secOrdAttr    <- "secOrd"              # Attribute storing a boolean whether the variational matrix has second order derivatives (TRUE) or not (FALSE).
etaAttr       <- "eta"                 # Attribute where eta parameter values of variational matrix are stored. The attribute contains a named numeric vector with the eta values. The names are the parameter names.
epsAttr       <- "eps"                 # Attribute where eps parameter values of variational matrix are stored. The attribute contains a named numeric vector with the eps values. The names are the parameter names.
parNormAttr   <- "parameterNormalized" # Attribute where parameter normalization status of variational matrix is stored. The attribute contains a named boolean vector with TRUE = normalized, FALSE = not. The names are the parameter names.
outNormAttr   <- "outputNormalized"    # Attribute where output normalization status of variational matrix is stored. The attribute contains a single boolean with TRUE = normalized, FALSE = not.

varTypeFim    <- "FIM"                 # Value of "type" attribute of FIM
omegaAttr     <- "omega"               # Attribute where omega parameter values of FIM are stored. The attribute contains a named numeric vector with the omega values. The names are the parameter names.
sigmaAttr     <- "sigma"               # Attribute where sigma parameter values of FIM are stored. The attribute contains a named numeric vector with the sigma values. The names are the parameter names.
thetaNormAttr <- "thetaNormalized"     # Attribute where theta parameter normalization status of FIM is stored. The attribute contains a named boolean vector with TRUE = normalized, FALSE = not. The names are the parameter names.
omegaNormAttr <- "omegaNormalized"     # Attribute where theta parameter normalization status of FIM is stored. The attribute contains a named boolean vector with TRUE = normalized, FALSE = not. The names are the parameter names.
sigmaNormAttr <- "sigmaNormalized"     # Attribute where theta parameter normalization status of FIM is stored. The attribute contains a named boolean vector with TRUE = normalized, FALSE = not. The names are the parameter names.

#------------------------------------------ processWarns ------------------------------------------
#' Process warnings
#'
#' Process a vector of warnings, by printing them to stdout.
#'
#' @param warns    A vector of warnings, of the form <caller>: <message>
#'
#' @return \code{NULL}
#'
#' @author Martijn van Noort
#' 
#' @noRd 
processWarns <- function(warns) {
  if (length(warns) > 0) cat(paste0("Warning in ", warns), sep = "")
}


#------------------------------------------ processErrors ------------------------------------------
#' Process errors
#'
#' Process a vector of errors, by printing them to stdout.
#'
#' @param errors    A vector of errors, of the form <caller>: <message>
#'
#' @return \code{NULL}
#'
#' @author Martijn van Noort
#' 
#' @noRd 
processErrors <- function(errors) {
  if (length(errors) > 0) cat(paste0("Error in ", errors), sep = "")
}


#--------------------------------------------------------------------------------------------------------
#------------------------------------------ Exported functions ------------------------------------------
#--------------------------------------------------------------------------------------------------------

#------------------------------------------ [.keepattr ------------------------------------------
#' Extract or replace parts of an object, keeping attributes
#'
#' This is a special version for objects (in particular data frames or matrices) with the 'keepattr' class,
#' that will keep attributes on extraction.
#' If the object \code{x} is a variational matrix or FIM, they will be kept only if the result is still a valid
#' variational matrix or FIM, respectively.
#' 
#' @param x   object from which to extract element(s) or in which to replace element(s).
#' @param ... other arguments passed to a parent generic. Typically used to specify a subset.
#'
#' @return The specified subset of \code{x}. If \code{x} is a variational matrix or FIM, and the result
#'   of subsetting is not a valid variational matrix or FIM, then the attributes will be lost.
#'   This can happen for example when the FIM is no longer a square matrix.
#'   They will also be lost if the class changes.
#'   This can happen for example when a single element is extracted.
#'
#' @export
#' 
#' @family retrieval functions
#'
#' @author Martijn van Noort
`[.keepattr` <- function(x, ...) {
  attrs <- attributes(x)
  attrs <- attrs[!names(attrs) %in% c("names", "dim", "dimnames", "row.names", "class", "tsp")]  # only keep non-special attributes and "comment"
  out <- x
  cls <- class(out) <- setdiff(class(x), "keepattr")  # remove keepattr so that data.frame subsetter will be used on next line
  out <- out[...]
  x <- out   # below, out will get attributes, x will be the version without
  if (identical(class(out), cls)) {
    # restore keepattr and other attributes, but only if other classes remained the same:
    for (nam in names(attrs)) attr(out, nam) <- attrs[[nam]]                 # restore attributes
    class(out) <- c("keepattr", class(out))
  }
  if(typeAttr %in% names(attributes(out)) && attr(out, typeAttr) == varTypeVar) {
    # Variational matrix. Check validity:
    valid <- isValidVariations(df = out, nmdf = "input")
    if (!is.null(valid)) {
      return(x)
    }
  }
  if(typeAttr %in% names(attributes(out)) && attr(out, typeAttr) == varTypeFim) {
    # FIM. Check validity:
    valid <- isValidFim(df = out, nmdf = "input")
    if (!is.null(valid)) {
      return(x)
    }
  }
  return(out)
}


#------------------------------------------ print.keepattr ------------------------------------------
#' Print values
#'
#' This is a special version for objects (in particular data frames or matrices) with the 'keepattr' class.
#' It removes the 'keepattr' class and prints according to the remaining class specification.
#' If \code{x} is a variational matrix or FIM, it is printed without attributes. 
#' 
#' @param x   object to be printed.
#' @param ... other arguments passed to the print method.
#'
#' @return \code{x}
#'
#' @export
#' 
#' @family plotting and printing
#'
#' @author Martijn van Noort
print.keepattr <- function(x, ...) {
  out <- x
  class(out) <- setdiff(class(x), "keepattr")
  if(typeAttr %in% names(attributes(out)) && attr(out, typeAttr) == varTypeVar) {
    print.data.frame(out)
  } else if(typeAttr %in% names(attributes(out)) && attr(out, typeAttr) == varTypeFim) {
    print.table(out)
  } else {
    print(out)
  }
  return(invisible(x))
}

