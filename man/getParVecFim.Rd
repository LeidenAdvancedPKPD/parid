% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FisherInfo.r
\name{getParVecFim}
\alias{getParVecFim}
\title{Get the column parameter values of a FIM as a single vector}
\usage{
getParVecFim(fim)
}
\arguments{
\item{fim}{The FIM.}
}
\value{
A named numeric vector containing the theta's, omega's and sigma's used to
create the FIM, and present as columns of the FIM, in the order of appearance in the FIM.
The names are the names of the parameters, where off-diagonal elements of omega and sigma are
named "(<par1>,<par2>)".
The function displays an error and returns \code{NULL} if \code{fim} is not a valid FIM.
}
\description{
Gets the values of the parameters corresponding to the columns of the FIM.
This is a subset of \code{getAllParsFim} corresponding to the columns present in the FIM.
}
\seealso{
Other retrieval functions: 
\code{\link{[.keepattr}()},
\code{\link{getAllParsFim}()},
\code{\link{getAllParsVariations}()},
\code{\link{getParVecVariations}()},
\code{\link{getParsFim}()},
\code{\link{getParsVariations}()},
\code{\link{subFim}()}
}
\author{
Martijn van Noort
}
\concept{retrieval functions}
