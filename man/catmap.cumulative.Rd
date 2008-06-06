\name{catmap.cumulative} 
\alias{catmap.cumulative}
\title{Cumulative Meta-Analyses and Plots using either Fixed- or Random-Effects} 
\author{Kristin K. Nicodemus, \email{kristin.nicodemus@well.ox.ac.uk}}
\description{ 
\code{catmap.cumulative} conducts cumulative meta-analyses and creates plots of the ORs and confidence intervals using either fixed- or random-effects analyses, and saves text files and plot files to the current working directory (use getwd() to obtain the current working directory). The plots are not created in the R Graphics device.
} 
\usage{catmap.cumulative(catmapobject, fe.cumulative, re.cumulative, fe.cumplot, 
re.cumplot) } 
\arguments{
	\item{catmapobject}{The catmap object created by a previous call to catmap}   
	\item{fe.cumulative}{Logical.  Should a cumulative meta-analysis be performed using fixed-effects estimates?  catmap assumes the order in which the studies are listed is the chronological ordering.  Automatic output result file is saved with the default name of \bold{dataset.fixed.effects.cumulative
.txt}, where dataset is the name of the file given as the first argument to catmap.  Note that repeated runs of the same input file will be appended to the default output file.}
	\item{re.cumulative}{Logical.  Should a cumulative meta-analysis be performed using random-effects estimates?  catmap assumes the order in which the studies are listed is the chronological ordering.  Automatic output result file is saved with the default name of \bold{dataset.random.effects.cumulative .txt}, where dataset is the name of the file given as the first argument to catmap.  Note that repeated runs of the same input file will be appended to the default output file.  Also note that random-effects estimates are undefined for a single study; calculations begin with the first 2 studies and then adds studies to them in a cumulative fashion.  The OR and CI for the first study may be found using the fixed effects estimates.}
	\item{fe.cumplot}{Logical.  Should a .pdf plot of the ORs and CIs from the cumulative meta-analysis using fixed-effects be output?  Can be TRUE only if fe.cumulative =TRUE.  Output plot file is saved with the default name of \bold{dataset.fixed. effects.cumulative.plot.pdf} where dataset is the name of the file given as the first argument to catmap.}
	\item{re.cumplot}{Logical.  Should a .pdf plot of the ORs and CIs from the cumulative meta-analysis using random-effects be output?  Can be TRUE only if re.cumulative= TRUE.  Output plot file is saved with the default name of \bold{dataset.random. effects
.cumulative.plot.pdf} where dataset is the name of the file given as the first argument to catmap.}
}
\details{
\code{catmap.cumulative} conducts cumulative meta-analyses and creates .pdf files of plots of the ORs and CIs using either the fixed-effect or the random-effect estimates.  \bold{NOTE: The studies should be listed in chronological order in the input file.  \code{catmap.cumulative} does not re-order studies by publication year.}
}
\seealso{
  \code{\link{catmap}},
  \code{\link{catmap.forest}},
  \code{\link{catmap.sense}},
  \code{\link{catmap.funnel}}.
}   
\examples{
\dontrun{
data(catmapdata) 	
catmapobject1<-catmap(catmapdata, 0.95, TRUE)
catmap.cumulative(catmapobject1, TRUE, TRUE, TRUE, TRUE)}}
\keyword{methods}
