% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VariationalEq.r
\name{getParVecVariations}
\alias{getParVecVariations}
\title{Get column parameter values of a variational matrix as a single vector}
\usage{
getParVecVariations(df)
}
\arguments{
\item{df}{A data frame containing a variational matrix, optionally with second order derivatives.}
}
\value{
A named numeric vector containing the theta's, eta's and eps's used to
create the variational matrix and present as columns of the matrix, in the order of appearance in the matrix.
The names are the names of the parameters.
The function displays an error and returns \code{NULL} if \code{df} is not a valid variational matrix.
}
\description{
Gets the parameter values of the columns of a variational matrix that were used to create a given variational matrix,
as a single vector.
This is the subset of \code{getAllParsVariations} corresponding to the columns present in the variational matrix.
The output is the same as that of \code{getParVariations} except it is formatted as a single vector.
}
\seealso{
Other retrieval functions: 
\code{\link{[.keepattr}()},
\code{\link{getAllParsFim}()},
\code{\link{getAllParsVariations}()},
\code{\link{getParVecFim}()},
\code{\link{getParsFim}()},
\code{\link{getParsVariations}()},
\code{\link{subFim}()}
}
\author{
Martijn van Noort
}
\concept{retrieval functions}
