\name{catmap} 
\alias{catmap}
\docType{package}
\title{catmap: Case-control And TDT Meta-Analysis Package} 
\author{Kristin K. Nicodemus, \email{nicodemusk@mail.nih.gov} or \email{knicodem@jhsph.edu}}
\description{ 
\code{catmap} is an R package that conducts fixed-effects (inverse variance) and random-effects (DerSimonian and Laird, 1986) meta-analyses of case-control or family-based (TDT) genetic data; in addition, it performs meta-analyses combining these two types of study designs.  The fixed-effects model was first described by Kazeem and Farrell (2005); the random-effects model is described in Nicodemus (submitted).} 
\usage{catmap(dataset, 
	ci, 
	fe.forest, 
	re.forest, 
	printout, 
	fe.sense, 
	fe.senseplot, 
	re.sense, 
	re.senseplot, 
	fe.cumulative, 
	fe.cumplot, 
	re.cumulative, 
	re.cumplot, 
	funnel) } 
\arguments{
	\item{dataset}{A text file containing a header with the following column names: \bold{name, study, t, nt, caserisk, controlrisk, casenotrisk, controlnotrisk} in tab-delimited format.  Note that the header must be exactly as specified and that all cells in the table must have an entry, even if the entry is 0 or missing (denoted NA).  See for example: data(catmapdata).)}   
	\item{ci}{The confidence level for confidence intervals; 0 < ci < 1}   
	\item{fe.forest}{Logical.  Should a postscript of the forest plot for the fixed-effects model be generated?  Plots are saved to the current working directory and have the default name of \bold{dataset.fixed.effects.plot.ps} where dataset is the name of the file given as the first argument to catmap.}
	\item{re.forest}{Logical.  Should a postscript of the forest plot for the random-effects model be generated?  Plots are saved to the current working directory and have the default name of \bold{dataset.random.effects.plot.ps} where dataset is the name of the file given as the first argument to catmap.}
	\item{printout}{Logical.  Should a text file of the fixed- and random-effects models, Q statistic, publication bias results be saved to the current working directory?  Output files are saved with the default name of \bold{dataset.output.txt} where dataset is the name of the file given as the first argument to catmap.}
	\item{fe.sense}{Logical.  Should a leave-one-out sensitivity analysis be performed using fixed-effects estimates?  Automatic output result files is saved with the default name of \bold{dataset.fixed.effects.sensitivity.txt}, where dataset is the name of the file given as the first argument to catmap.  Note that repeated runs of the same input file will be appended to the default output file names; only the output option overwrites existing files.}
	\item{fe.senseplot}{Logical.  Should a .pdf plot of the ORs and CIs from the sensitivity analysis using fixed-effects be output?  Can be TRUE only if fe.sense=TRUE.  Output plot file is saved with the default name of \bold{dataset.fixed.effects.sensitivity.plot.pdf} where dataset is the name of the file given as the first argument to catmap.}
	\item{re.sense}{Logical.  Should a leave-one-out sensitivity analysis be performed using random-effects estimates?  Automatic output result file is saved with the default name of \bold{dataset.random.effects.sensitivity.txt}, where dataset is the name of the file given as the first argument to catmap.  Note that repeated runs of the same input file will be appended to the default output file names; only the output option overwrites existing files.}
	\item{re.senseplot}{Logical.  Should a .pdf plot of the ORs and CIs from the sensitivity analysis using random-effects be output?  Can be TRUE only if re.sense=TRUE.  Output plot file is saved with the default name of \bold{dataset.random.effects.sensitivity.plot. pdf} where dataset is the name of the file given as the first argument to catmap.}
	\item{fe.cumulative}{Logical.  Should a cumulative meta-analysis be performed using fixed-effects estimates?  catmap assumes the order in which the studies are listed is the chronological ordering.  Automatic output result file is saved with the default name of \bold{dataset.fixed.effects.cumulative.txt}, respectively, where dataset is the name of the file given as the first argument to catmap.  Note that repeated runs of the same input file will be appended to the default output file names; only the output option overwrites existing files.}
	\item{fe.cumplot}{Logical.  Should a .pdf plot of the ORs and CIs from the cumulative meta-analysis using fixed-effects be output?  Can be TRUE only if fe.cumulative=TRUE.  Output plot file is saved with the default name of \bold{dataset.fixed.effects. cumulative.plot.pdf} where dataset is the name of the file given as the first argument to catmap.}
	\item{re.cumulative}{Logical.  Should a cumulative meta-analysis be performed using random-effects estimates?  catmap assumes the order in which the studies are listed is the chronological ordering.  Automatic output result file is saved with the default name of \bold{dataset.random.effects.cumulative.txt}, where dataset is the name of the file given as the first argument to catmap.  Note that repeated runs of the same input file will be appended to the default output file names; only the output option overwrites existing files.  Also note that random-effects estimates are undefined for a single study; calculations begin with the first 2 studies and adds studies to them in a cumulative fashion.  The OR and CI for the first study may be found using the fixed effects estimates.}
	\item{re.cumplot}{Logical.  Should a .pdf plot of the ORs and CIs from the cumulative meta-analysis using random-effects be output?  Can be TRUE only if re.cumulative= TRUE.  Output plot file is saved with the default name of \bold{dataset.random.effects.
cumulative.plot.pdf} where dataset is the name of the file given as the first argument to catmap.}
	\item{funnel}{Logical.  Should a funnel plot be produced?  Funnel plots plot the log standard error of the ORs and the ORs along with a solid line at 1.0 (no effect) and a dotted line at the overall pooled OR.  Used to assess publication bias.  Output plot file is saved with the default name of \bold{dataset.funnel.plot.pdf} where dataset is the name of the file given as the first argument to catmap.}
} 
\details{catmap is an R package that conducts fixed-effects (inverse variance) and random-effects (DerSimonian and Laird, 1986) meta-analyses of case-control or family-based (TDT) genetic data; in addition, it performs meta-analyses combining these two types of study designs.  The fixed-effects model was first described by Kazeem and Farrell (2005); the random-effects model is described in Nicodemus (submitted). Cumulative meta-analyses over time and leave-one-out sensitivity analyses may be performed using either fixed- or random-effects estimates or both estimates may be calculated; both produce a .txt file and an optional .pdf plot as output. A funnel plot graphic is implemented; however, no formal test of publication bias is available (see Ioannidis & Trikalinos, 2007).} 
\examples{
data(catmapdata) 	
catmap(catmapdata, 0.95, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE)
\dontrun{
catmap(catmapdata, 0.95, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE)}
}  
\keyword{methods}

