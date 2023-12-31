% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sensitivity.r
\name{calcSensitivityFromModel}
\alias{calcSensitivityFromModel}
\title{Sensitivity indices, parameter relations and least identifiable parameters}
\usage{
calcSensitivityFromModel(
  outputs,
  model,
  parms,
  init,
  outputPred,
  times,
  vars = names(parms),
  parNorm = FALSE,
  outNorm = FALSE,
  symbolic = TRUE,
  chkModel = TRUE,
  Mthresh = NULL,
  normalizeM = TRUE,
  scaleM = TRUE,
  ...
)
}
\arguments{
\item{outputs}{List of the indices and parameter relations to be determined.
List of strings, where each element is one of "N", "S", "A", "R", "M", "L".
They stand for null space dimension, skewing coefficient, skewing angle, parameter relations,
minimal parameter relations and least identifiable parameters, respectively.}

\item{model}{Function(t, y, p) of time, state variables and parameters, specifying the differential equations.
This function should return the numeric vector dy/dt.}

\item{parms}{Model parameter values, as named numeric vector.}

\item{init}{Function(p) of parameters specifying the initial state of the model.}

\item{outputPred}{Function(y, p) of state variables and parameters specifying the model prediction, i.e., the
output without residual.
This function should return a numeric vector.}

\item{times}{Numeric vector of times where variations are to be evaluated.}

\item{vars}{Vector of names of variables q to be included in sensitivity analysis.
Should be a subset of \code{names(parms)}.
By default equal to \code{names(parms)}.
Variations d(outputPred)/dq are calculated only for these parameters.
If specified as \code{NULL}, it is set to \code{names(parms)}.}

\item{parNorm}{Vector (not necessarily named) of booleans of the same dimension as \code{parms}.
The value \code{TRUE} means the corresponding individual parameter should be normalized in the FIM,
\code{FALSE} means not.
If the vector is named, the names should match \code{vars}.
Instead of a vector, may also provide a single boolean, where \code{TRUE} stands for a vector of
all \code{TRUE}, and \code{FALSE} (default) for a vector of all \code{FALSE}.}

\item{outNorm}{\code{TRUE} if variations should be output-normalized, \code{FALSE} (default) if not.}

\item{symbolic}{\code{TRUE} (default) if derivatives are to be computed symbolically, \code{FALSE} if numerically.
See \code{\link{calcVariations}} for notes on using symbolic derivation.}

\item{chkModel}{\code{TRUE} (default) if it has to be checked whether model components (\code{model}, \code{p}, \code{init}, \code{output})
are formatted correctly for symbolic manipulation, \code{FALSE} if not.}

\item{Mthresh}{Threshold for minimal parameter relations.
Parameter relations smaller than \code{Mthresh/sqrt(k)} (before scaling, see \code{scaleM}) are returned,
where k is the number of parameters.
The factor sqrt(k) makes Mthresh independent of dimension.
If \code{NULL}, only the most minimal relation is returned.}

\item{normalizeM}{\code{TRUE} (default) if the minimal parameter relations should be normalized, \code{FALSE} if not.
This is not the same as parameter- or output-normalization of the sensitivity matrix.
See \code{\link{calcSensitivityFromMatrix}} for more information.}

\item{scaleM}{\code{TRUE} (default) if norms of the minimal parameter relations should be scaled to the range \link{0, 1},
\code{FALSE} if not.
If \code{TRUE}, the norms are divided by \code{sqrt(k)}, where k is the number of parameters.
If \code{FALSE} they are unchanged.}

\item{...}{Named arguments to be passed to \code{\link[deSolve]{lsoda}}.
Can be used for example to pass events or error tolerances.}
}
\value{
A list with a named element for each element of 'outputs'. Elements for "N", "S", "A" are numeric values,
while elements "R", "M" contain row-named matrices with the parameter relations defined by its columns,
and "L" contains a named list with one element for each parameter, in increasing order of identifiability or distance,
of the form "parameter name" = "distance of this parameter to the span of the others" in observation space.
Distance (near) zero indicates non-identifiability of this parameter.
Element M has an attribute "norm" specifying the norms of the minimal parameter relations, where near zero values
correspond to non-identifiability, and large values to identifiability.
The relations are provided in order of increasing norm (if there is more than one).
It also has attributes "normalized" and "scaled" containing the values of the input parameters \code{normalizeM}
and \code{scaleM}, respectively.

If the variational matrix could not be computed or \code{vars} is not defined correctly,
then the return value is \code{NULL} and an error message is shown.
If there is a single parameter in \code{vars}, then elements S, A and L are \code{NA}.
}
\description{
Calculates sensitivity indices, parameter relations and/or the least identifiable parameters for a given model.
It first computes the variational matrix for the model and then applies \code{\link{calcSensitivityFromMatrix}}.
}
\note{
Output-normalization divides the vector dy_i/dp_j by y_i, and is recommended if the output is
considered to be a multiplicative rather than an additive variable.
Parameter-normalization multiplies the vector dy_i/dp_j by p_j, and is recommended if the input is
considered to be a multiplicative rather than an additive variable.
(A multiplicative variable is a variable where differences between values are considered on the log scale.
So if e.g. values 0.001 and 0.01 are considered different because they differ by a factor of 10, then the
variable is multiplicative.
If they are considered approximately equal because they are both close to 0,
then the variable is additive.
Most PK inputs and outputs are multiplicative, likewise for many PD ones.)

Output-normalization is equivalent to considering the sensitivity matrix dlog(y)/dp rather than dy/dp.
Parameter-normalization is equivalent to considering the sensitivity matrix dy/dlog(p) rather than dy/dp.

It is recommended to use symbolic derivation (\code{symbolic==TRUE}) rather than numerical derivation,
as it is more accurate, especially for derivatives that are formally zero.
See \code{\link{calcVariations}} for details.
}
\seealso{
\code{\link{calcSensitivityFromMatrix}} for details on the calculation of the outputs.

Other calculation functions: 
\code{\link{calcAliasingScoresFromMatrix}()},
\code{\link{calcAliasingScoresFromModel}()},
\code{\link{calcFimFromMatrix}()},
\code{\link{calcFimFromModel}()},
\code{\link{calcFirstVariations}()},
\code{\link{calcSensitivityFromMatrix}()},
\code{\link{calcVariationsFim}()},
\code{\link{calcVariations}()},
\code{\link{fimIdent}()}
}
\author{
Martijn van Noort
}
\concept{calculation functions}
