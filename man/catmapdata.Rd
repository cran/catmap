\name{catmapdata}
\alias{catmapdata}
\docType{data}
\title{An example data frame for use with catmap}
\description{
A file for use with catmap containing simulated data.
}
\usage{data(catmapdata)}
\format{
  A data frame with 5 observations on the following 8 variables.
  \describe{
    \item{\code{name}}{a factor with study name and optionally year of publication.  NOTE: if year of publication is included there must be no space between study name and year.  A comma or underscore work nicely. \code{Abrams,2001} \code{Peter,2002} \code{Todd,2003} \code{Wei,2007} \code{Yu,2007}}
    \item{\code{study}}{a numeric vector containing 1 if study is TDT and 2 if case-control}
    \item{\code{t}}{a numeric vector containing counts of alleles transmitted in TDT study}
    \item{\code{nt}}{a numeric vector containing counts of alleles not transmitted in TDT study}
    \item{\code{caserisk}}{a numeric vector containing counts of risk alleles in cases}
    \item{\code{controlrisk}}{a numeric vector containing counts of risk alleles in controls}
    \item{\code{casenotrisk}}{a numeric vector containing counts of non-risk alleles in cases}
    \item{\code{controlnotrisk}}{a numeric vector containing counts of non-risk alleles in controls}
  }
}
\details{
The header must be part of the file and either 0 or NA must be included for data not relevant for the particular study.  For example, using a TDT study the caserisk, controlrisk, casenotrisk and controlnotrisk must have values of either 0 or NA.
}
\examples{
data(catmapdata)
}
\keyword{datasets}
