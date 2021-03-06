% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{normFactors}
\alias{normFactors}
\alias{normFactors,srnadiffExp-method}
\alias{normFactors<-}
\alias{normFactors<-,srnadiffExp,numeric-method}
\title{Accessors for the 'normFactors' slot of an srnadiffExp object}
\usage{
\S4method{normFactors}{srnadiffExp}(object)

\S4method{normFactors}{srnadiffExp,numeric}(object) <- value
}
\arguments{
\item{object}{An \code{srnadiffExp} object.}

\item{value}{A numeric vector, one size factor for each sample
in the coverage data.}
}
\value{
The normalization factors, in a list.
}
\description{
The \code{normFactors} slot holds the normalization factors as a
named vector with sample names.
}
\details{
The \code{normFactors} vector assigns to each sample coverage a value,
the normalization factor, such that count values in each sample coverage
can be brought to a common scale by dividing by the corresponding
normalization factor. This step is also called normalization, its purpose
is to render coverages (counts) from different samples, which may have
been sequenced to different depths, comparable. Normalization factors
are estimated using the \emph{median ratio method} described by
Equation 5 in Anders and Huber (2010). Alternative normalization factor
estimators can also be supplied using the assignment function
\code{sizeFactors<-}.
}
\examples{
srnaExp <- srnadiffExample()
normFactors(srnaExp)

}
\references{
Simon Anders, and Wolfgang Huber (2010). Differential expression analysis
for sequence count data. \emph{Genome Biology}, 11:106.
}
