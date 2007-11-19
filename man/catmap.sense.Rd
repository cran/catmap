\name{catmap.sense} 
\alias{catmap.sense}
\title{Leave-One-Out Sensitivity Analyses and Plots using either Fixed- or Random-Effects Estimates} 
\author{Kristin K. Nicodemus, \email{kristin.nicodemus@well.ox.ac.uk}}
\description{ 
\code{catmap.sense} conducts leave-one-out sensitivity analyses and creates plots of the ORs and confidence intervals using either fixed- or random-effects analyses, which are saved to the current working directory (use getwd() to view) but does not create the plots in the R Graphics device.
} 
\usage{catmap.sense(catmapobject, fe.sense, re.sense, fe.senseplot, re.senseplot) } 
\arguments{
	\item{catmapobject}{The catmap object created by a previous call to catmap}   
	\item{fe.sense}{Logical.  Should a leave-one-out sensitivity analysis be performed using fixed-effects estimates?  Automatic output result files are saved with the default name of \bold{dataset.fixed.effects.sensitivity.txt}, where dataset is the name of the file given as the first argument to catmap.  Note that repeated runs of the same input file will be appended to the default output file names.}
	\item{re.sense}{Logical.  Should a leave-one-out sensitivity analysis be performed using random-effects estimates?  Automatic output result files are saved with the default name of \bold{dataset.random.effects.sensitivity.txt}, where dataset is the name of the file given as the first argument to catmap.  Note that repeated runs of the same input file will be appended to the default output file names.}
	\item{fe.senseplot}{Logical.  Should a .pdf plot of the ORs and CIs from the sensitivity analysis using fixed-effects be output?  Can be TRUE only if fe.sense=TRUE.  Output plot file is saved with the default name of \bold{dataset.fixed.effects.sensitivity.plot.pdf} where dataset is the name of the file given as the first argument to catmap.}
	\item{re.senseplot}{Logical.  Should a .pdf plot of the ORs and CIs from the sensitivity analysis using random-effects be output?  Can be TRUE only if re.sense=TRUE.  Output plot file is saved with the default name of \bold{dataset.random.effects.sensitivity.plot. pdf} where dataset is the name of the file given as the first argument to catmap.}
}
\details{
\code{catmap.sense} conducts leave-one-out sensitivity analyses and creates .pdf files of plots of the ORs and CIs using either the fixed-effect or the random-effect estimates. Plots are not created in the R Graphics device window, but are instead saved to a .pdf file in the current working directory. Likewise the .txt files of results will be saved to the current working directory. To find the current working directory, use getwd()
}
\seealso{
  \code{\link{catmap}},
  \code{\link{catmap.forest}},
  \code{\link{catmap.cumulative}},
  \code{\link{catmap.funnel}}.
}  
\examples{
\dontrun{
data(catmapdata) 	
catmapobject1<-catmap(catmapdata, 0.95, TRUE)
catmap.sense(catmapobject1, TRUE, TRUE, TRUE, TRUE)}}
\keyword{methods}