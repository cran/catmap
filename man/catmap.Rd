\name{catmap} 
\alias{catmap}
\docType{package}
\title{catmap: Case-control And TDT Meta-Analysis Package} 
\author{Kristin K. Nicodemus, \email{kristin.nicodemus@well.ox.ac.uk}}
\description{ 
\code{catmap} is an R package that conducts fixed-effects (inverse variance) and random-effects (DerSimonian and Laird, 1986) meta-analyses of case-control or family-based (TDT) genetic data; in addition, it performs meta-analyses combining these two types of study designs.  The fixed-effects model was first described by Kazeem and Farrall (2005); the random-effects model is described in Nicodemus (submitted) and saves a text file to the current working directory of all results printed to screen (use getwd() to find cwd).
} 
\usage{catmapobject<-catmap(dataset, ci, printout) } 
\arguments{
	\item{dataset}{A text file containing a header with the following column names: \bold{name, study, t, nt, caserisk, controlrisk, casenotrisk, controlnotrisk} in tab-delimited format.  Note that the header must be exactly as specified and that all cells in the table must have an entry, even if the entry is 0 or missing (NA).  See for example: data(catmapdata). \bold{The dataset argument to catmap should be either the example data or a file containing the data for catmap, not an R object.}
}   
	\item{ci}{The confidence level for confidence intervals; 0 < ci < 1.  The default is 0.95}   
	\item{printout}{Logical.  Should a text file of the fixed- and random-effects models and Q statistic results be saved to the current working directory?  Output files are saved with the default name of \bold{dataset.output.txt} where dataset is the name of the file given as the first argument to catmap.  Default = TRUE.}
}
\details{catmap is an R package that conducts fixed-effects (inverse variance) and random-effects (DerSimonian and Laird, 1986) meta-analyses of case-control or family-based (TDT) genetic data; in addition, it performs meta-analyses combining these two types of study designs.  The fixed-effects model was first described by Kazeem and Farrall (2005); the random-effects model is described in Nicodemus (submitted). Cumulative meta-analyses over time and leave-one-out sensitivity analyses may be performed using either fixed- or random-effects estimates or both estimates may be calculated; both produce a .txt file and an optional .pdf plot as output. A funnel plot graphic is implemented; however, no formal test of publication bias is available (see Ioannidis & Trikalinos, 2007). If users request a file to be created containing the results the file will be saved to the current working directory, which users can find by using >getwd(). \bold{Note that a catmap object must be created on the first call to catmap. If you do not create a catmap object you will not be able to use any of the other functions AND you will get a printout of the entire contents to screen} 
} 
\seealso{
  \code{\link{catmap.forest}},
  \code{\link{catmap.sense}},
  \code{\link{catmap.cumulative}},
  \code{\link{catmap.funnel}}.
}  
\examples{
\dontrun{
data(catmapdata) 	
catmap(catmapdata, 0.95, TRUE)
catmapobject1<-catmap(catmapdata, 0.95, TRUE)}}
\keyword{methods}

