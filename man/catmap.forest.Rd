\name{catmap.forest} 
\alias{catmap.forest}
\title{Forest Plots using either Fixed- or Random-Effects Pooled ORs and CIs} 
\author{Kristin K. Nicodemus, \email{kristin.nicodemus@well.ox.ac.uk}}
\description{ 
\code{catmap.forest} creates forest plots of the individual study ORs and CIs and the fixed or random effects pooled OR and CI and saves them to the current working directory (use getwd() to view cwd).  The plots are not created in the R Graphics device.
} 
\usage{catmap.forest(catmapobject, fe.forest, re.forest) } 
\arguments{
	\item{catmapobject}{The catmap object created by a previous call to catmap}   
	\item{fe.forest}{Logical.  Should a forest plot be created using the fixed-effects estimates?  Plots are saved with the default name of \bold{dataset.fixed.effects.forest.pdf}, where dataset is the name of the file given as the first argument to catmap.}
	\item{re.forest}{Logical.  Should a forest plot be created using the random-effects estimates?  Plots are saved with the default name of \bold{dataset.random.effects.forest.pdf} where dataset is the name of the file given as the first argument to catmap.}
}
\details{
\code{catmap.forest} creates forest plots of individual study ORs and CIs plus the pooled estimate of the fixed- or random-effects pooled OR and CI. Plots are not created in the R Graphics device window, but are instead saved to a .pdf file in the current working directory, which can be found using getwd()
} 
\seealso{
  \code{\link{catmap}},
  \code{\link{catmap.cumulative}},
  \code{\link{catmap.sense}},
  \code{\link{catmap.funnel}}.
}  
\examples{
\dontrun{
data(catmapdata) 	
catmapobject1<-catmap(catmapdata, 0.95, TRUE)
catmap.cumulative(catmapobject1, TRUE, TRUE, TRUE, TRUE)}}
\keyword{methods}
