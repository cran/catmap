\name{catmap.funnel} 
\alias{catmap.funnel}
\title{Funnel Plots for catmap} 
\author{Kristin K. Nicodemus, \email{nicodemusk@mail.nih.gov} or \email{knicodem@jhsph.edu}}
\description{ 
\code{catmap.funnel} creates a funnel plot of the individual ORs against the se.
} 
\usage{catmap.funnel(catmapobject, funnel) } 
\arguments{
	\item{catmapobject}{The catmap object created by a previous call to catmap}   
	\item{funnel}{Logical.  Should a funnel plot be produced?  Funnel plots plot standard error of the log ORs against the ORs along with a solid line at 1.0 and a dotted line at the overall pooled OR.  Used to assess publication bias.  Output plot file is saved with the default name of \bold{dataset.funnel.plot.pdf} where dataset is the name of the file given as the first argument to catmap.}
} 
\details{\code{catmap.funnel} creates a .pdf file of the funnel plot.
} 
\seealso{
  \code{\link{catmap}},
  \code{\link{catmap.sense}},
  \code{\link{catmap.cumulative}},
  \code{\link{catmap.forest}}.
}  
\examples{
\dontrun{
data(catmapdata) 	
catmapobject1<-catmap(catmapdata, 0.95, TRUE)
catmap.funnel(catmapobject1, TRUE)}}
\keyword{methods}
