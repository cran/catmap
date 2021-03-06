% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global.R
\docType{data}
\name{catmapdata}
\alias{catmapdata}
\title{Example \code{catmap} Data}
\format{A \code{data.frame} with 5 observations and 8 variables.
\itemize{
\item{name:}{ a factor with study name and optionally year of publication.
 NOTE: if year of publication is included there must be no space between
 study name and year. A comma or underscore works nicely (e.g.,
 \code{Abrams,2001} \code{Peter,2002} \code{Todd,2003}
 \code{Wei,2007} \code{Yu,2007})
}\item{study:}{ a numeric vector containing 1 if the study is TDT and 2 if the study is case-control
}\item{t:}{ a numeric vector containing counts of alleles transmitted in the TDT study
}\item{nt:}{ a numeric vector containing counts of alleles not transmitted in the TDT study
}\item{caserisk:}{ a numeric vector containing counts of risk alleles in cases
}\item{controlrisk:}{ a numeric vector containing counts of risk alleles in controls
}\item{casenotrisk:}{ a numeric vector containing counts of non-risk alleles in cases
}\item{controlnotrisk:}{ a numeric vector containing counts of non-risk alleles in controls
 }}}
\usage{
catmapdata
}
\description{
An example data set for use with \code{catmap}. All input data should
 have the header as part of the file and either 0 or NA values for entries
 not relevant to that particular study design. For example, TDT studies
 should have the caserisk, controlrisk, casenotrisk and controlnotrisk
 values set to either 0 or NA.
}
\seealso{
\code{\link{catmap}}, \code{\link{catmap.forest}},
 \code{\link{catmap.sense}}, \code{\link{catmap.cumulative}},
 \code{\link{catmap.funnel}}
}
\author{
Algorithm designed and implemented by Kristin K. Nicodemus.
 Code modified and updated by Thom Quinn.
}
\keyword{datasets}
