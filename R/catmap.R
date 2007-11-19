catmap<-function(dataset, ci, printout){
options(warn=-1)
data(catmapdata)

#make defaults for catmap of 0.95 for CI and printout = TRUE
if(missing(ci)){
ci<-0.95
}
if(missing(printout)){
printout<-TRUE
}

#read in data from function
if(dataset!=catmapdata){
a1<-read.table(dataset, header=T)
}
if(dataset==catmapdata){
a1<-catmapdata
}

#split data into tdt and case-control studies
tdt<-a1[a1$study==1,]
cc<-a1[a1$study==2,]

#fixed-effects estimates
#calculate tdt-specific log ORs, variances and weights

tdtlogOR<-log(tdt$t/tdt$nt)
tdtvar<-((1/tdt$t)+(1/tdt$nt))
tdtweight<-(1/tdtvar)
tdtse<-sqrt(tdtvar)

#calculate case-control specific log ORs, variances and weights

cclogOR<-log((cc$caserisk*cc$controlnotrisk)/(cc$casenotrisk*cc$controlrisk))
ccvar<-((1/cc$caserisk)+(1/cc$controlrisk)+(1/cc$casenotrisk)+(1/cc$controlnotrisk))
ccse<-sqrt(ccvar)
ccweight<-(1/ccvar)

#calculate combined log OR, variance, confidence interval and p-value

a1$lev<-as.factor(a1$study)

if(nlevels(a1$lev)==1 & a1$study[1]==1){
weight<-tdtweight
logOR<-tdtlogOR
seLogOR<-tdtse
comvarlogOR<-tdtvar
}
if(nlevels(a1$lev)==1 & a1$study[1]==2){
weight<-ccweight
logOR<-cclogOR
seLogOR<-ccse
comvarlogOR<-ccvar
}
if(nlevels(a1$lev)==2){
studyorder<-cbind(1:nrow(a1), a1$study)
studyorder1<-studyorder[order(studyorder[,2]),]
weight1<-c(tdtweight, ccweight)
logOR1<-c(tdtlogOR, cclogOR)
seLogOR1<-c(tdtse, ccse)
comvarlogOR1<-c(tdtvar, ccvar)
studyorder2<-cbind(studyorder1, weight1, logOR1, seLogOR1, comvarlogOR1)
studyorder3<-studyorder2[order(studyorder2[,1]),]
weight<-studyorder3[,3]
logOR<-studyorder3[,4]
seLogOR<-studyorder3[,5]
comvarlogOR<-studyorder3[,6]
}

combinedLogOR<-((sum(weight*logOR))/sum(weight))
combinedOR<-exp(combinedLogOR)
combinedSeLogOR<-(sqrt(1/sum(weight)))
combinedVarLogOR<-(1/sum(weight))
combinedChisq<-(((combinedLogOR-0)^2)/combinedVarLogOR)
combinedValue<-pchisq(combinedChisq, df=1)
combinedPvalue<-(1-combinedValue)

#get qnorm values
alpha<-(1-((1-ci)/2))
quantilenorm<-qnorm(alpha, 0, 1)

lbci<-exp(combinedLogOR-(quantilenorm*combinedSeLogOR))
ubci<-exp(combinedLogOR+(quantilenorm*combinedSeLogOR))
combinedCI<-c(lbci, ubci)
SeLogOR<-sqrt(comvarlogOR)
lbci.fe<-exp(logOR-(quantilenorm*SeLogOR))
ubci.fe<-exp(logOR+(quantilenorm*SeLogOR))

#calculate heterogeneity
het.df<-(nrow(a1)-1)
chisqHet<-(sum(weight*(((logOR-combinedLogOR)^2))))
combinedHetValue<-pchisq(chisqHet, df=het.df)
heterogeneityPvalue<-(1-combinedHetValue)

#DerSimonian and Laird random-effects estimates

tau2<-((chisqHet-het.df)/(sum(weight)-(sum(weight^2)/(sum(weight)))))

#print results - print to screen

if(tau2 <=0){
table.header<-c("Inverse Variance Fixed-Effects OR", "Inverse Variance Fixed-Effects Lower Bound CI", "Inverse Variance Fixed-Effects Upper Bound CI", "Inverse Variance Fixed-Effects Chi-Square", "Inverse Variance Fixed-Effects p-value", "Q Statistic (Heterogeneity) Chi-Square", "Q Statistic (Heterogeneity) p-value")
table.fill<-c(combinedOR, combinedCI, combinedChisq, combinedPvalue, chisqHet, heterogeneityPvalue)
results<-rbind(table.header, round(table.fill, digits=5))
cat("NOTICE: tau2 is less than or equal to 0;\n no random effects estimates will be calculated\n Pooled Estimates\n")
cat(results, sep="\n")
cat("\n")
}

if (tau2 > 0){
weight.dsl<-(1/(comvarlogOR+tau2))
logOR.dsl<-((sum(weight.dsl*logOR))/(sum(weight.dsl)))
OR.dsl<-exp(logOR.dsl)
seLogOR.dsl<-(1/(sqrt(sum(weight.dsl))))
varLogOR.dsl<-(1/sum(weight.dsl))
lbci.dsl<-exp(logOR.dsl-(quantilenorm*seLogOR.dsl))
ubci.dsl<-exp(logOR.dsl+(quantilenorm*seLogOR.dsl))
ci.dsl<-c(lbci.dsl, ubci.dsl)
chisq.dsl<-(((logOR.dsl-0)^2)/varLogOR.dsl)
value.dsl<-pchisq(chisq.dsl, df=1)
pvalue.dsl<-(1-value.dsl)

#print results - print to screen

table.header<-c("Inverse Variance Fixed-Effects OR", "Inverse Variance Fixed-Effects Lower Bound CI", "Inverse Variance Fixed-Effects Upper Bound CI", "Inverse Variance Fixed-Effects Chi-Square", "Inverse Variance Fixed-Effects p-value", "Q Statistic (Heterogeneity) Chi-Square", "Q Statistic (Heterogeneity) p-value", "DerSimonian & Laird Random-Effects OR", "DerSimonian & Laird Random-Effects Lower Bound CI", "DerSimonian & Laird Random-Effects Chi-Square", "DerSimonian & Laird Random-Effects p-value")
table.fill<-c(combinedOR, combinedCI, combinedChisq, combinedPvalue, chisqHet, heterogeneityPvalue, OR.dsl, ci.dsl, pvalue.dsl)
results<-rbind(table.header, round(table.fill, digits=5))
cat("Pooled Estimates\n")
cat(results, sep="\n")
cat("\n")

}

studyname<-a1$name

#print results to file dataset.output.txt

if (printout==TRUE & dataset!=catmapdata){
sink(paste(dataset, "output.txt", sep="."))
cat("Pooled Estimates\n")
cat(results, sep="\n")
sink()
}

if (printout==TRUE & dataset==catmapdata){
sink(paste("catmapdata", "output.txt", sep="."))
cat("Pooled Estimates\n")
cat(results, sep="\n")
sink()
}

if(tau2 <=0){
return(comvarlogOR, combinedLogOR, combinedOR, combinedSeLogOR, weight, logOR, combinedVarLogOR, combinedChisq, combinedValue, combinedPvalue, lbci, ubci, combinedCI, SeLogOR, lbci.fe, ubci.fe, het.df, chisqHet, combinedHetValue, heterogeneityPvalue, tau2, studyname, a1, quantilenorm, ci, dataset) 
}
if(tau2 > 0){
return(comvarlogOR, combinedLogOR, combinedOR, combinedSeLogOR, weight, logOR, combinedVarLogOR, combinedChisq, combinedValue, combinedPvalue, lbci, ubci, combinedCI, SeLogOR, lbci.fe, ubci.fe, het.df, chisqHet, combinedHetValue, heterogeneityPvalue, tau2, studyname, a1, quantilenorm, ci, dataset, weight.dsl, logOR.dsl, OR.dsl, seLogOR.dsl, varLogOR.dsl, lbci.dsl, ubci.dsl, ci.dsl, chisq.dsl, value.dsl, pvalue.dsl)
}
}

#sensitivity analysis

catmap.sense<-function(catmapobject, fe.sense, re.sense, fe.senseplot, re.senseplot){

#fixed-effects sensitivity

if (fe.sense==TRUE){
sfplot<-matrix(0,nrow(catmapobject$a1),3)
for(f in 1:nrow(catmapobject$a1)){
sf.weight<-catmapobject$weight[-f]
sf.logOR<-catmapobject$logOR[-f]

sf.combinedLogOR<-((sum(sf.weight*sf.logOR))/sum(sf.weight))
sf.combinedOR<-exp(sf.combinedLogOR)
sf.combinedSeLogOR<-(sqrt(1/sum(sf.weight)))
sf.combinedVarLogOR<-(1/sum(sf.weight))
sf.combinedChisq<-(((sf.combinedLogOR-0)^2)/sf.combinedVarLogOR)
sf.combinedValue<-pchisq(sf.combinedChisq, df=1)
sf.combinedPvalue<-(1-sf.combinedValue)

#get qnorm values
alpha<-(1-((1-catmapobject$ci)/2))
quantilenorm<-qnorm(alpha, 0, 1)

sf.lbci<-exp(sf.combinedLogOR-(quantilenorm*sf.combinedSeLogOR))
sf.ubci<-exp(sf.combinedLogOR+(quantilenorm*sf.combinedSeLogOR))
sf.combinedCI<-c(sf.lbci, sf.ubci)
sf.SeLogOR<-sqrt(sf.combinedVarLogOR)
sf.lbci<-exp(sf.logOR-(quantilenorm*sf.SeLogOR))
sf.ubci<-exp(sf.logOR+(quantilenorm*sf.SeLogOR))

#calculate heterogeneity
sf.chisqHet<-(sum(sf.weight*(((sf.logOR-sf.combinedLogOR)^2))))
sf.df<-(nrow(catmapobject$a1)-2)
sf.combinedHetValue<-pchisq(sf.chisqHet, df=sf.df)
sf.heterogeneityPvalue<-(1-sf.combinedHetValue)

study.removed<-paste("Study Removed =", catmapobject$studyname[f], sep=" ")
sftable.header<-c("Inverse Variance Fixed-Effects OR", "Inverse Variance Fixed-Effects Lower Bound CI", "Inverse Variance Fixed-Effects Upper Bound CI", "Inverse Variance Fixed-Effects Chi-Square", "Inverse Variance Fixed-Effects p-value", "Q Statistic (Heterogeneity) Chi-Square", "Q Statistic (Heterogeneity) p-value")
sftable.fill<-c(sf.combinedOR, sf.combinedCI, sf.combinedChisq, sf.combinedPvalue, sf.chisqHet, sf.heterogeneityPvalue)
sf.results<-rbind(sftable.header, round(sftable.fill, digits=5))
sfvalues<-c(sf.combinedOR, sf.combinedCI)
sfplot[f,]<-sfvalues
cat("Fixed Effects Sensitivity Analysis\n")
cat(study.removed, sf.results, sep="\n")
cat("\n")

#print results to file dataset.fixed.effects.sensitivity.txt

if(catmapobject$dataset!=catmapdata){
sink(paste(catmapobject$dataset, "fixed.effects.sensitivity.txt", sep="."), append=TRUE)
cat("Fixed Effects Sensitivity Analysis\n")
cat(study.removed, sf.results, sep="\n")
cat("\n")
sink()
}

if(catmapobject$dataset==catmapdata){
sink(paste("catmapdata", "fixed.effects.sensitivity.txt", sep="."), append=TRUE)
cat("Fixed Effects Sensitivity Analysis\n")
cat(study.removed, sf.results, sep="\n")
cat("\n")
sink()
}
}
}

#random-effects sensitivity

if (re.sense==TRUE & catmapobject$tau2 <=0){
cat("NOTICE: tau2 is less than or equal to 0;\n no random effects estimates will be calculated.\n")
}

if (re.sense==TRUE & catmapobject$tau2 > 0){
srplot<-matrix(0,nrow(catmapobject$a1),3)
for(r in 1:nrow(catmapobject$a1)){
sr.weight<-catmapobject$weight[-r]
sr.logOR<-catmapobject$logOR[-r]
sr.comvarlogOR<-catmapobject$comvarlogOR[-r]

sr.combinedLogOR<-((sum(sr.weight*sr.logOR))/sum(sr.weight))
sr.combinedOR<-exp(sr.combinedLogOR)
sr.combinedSeLogOR<-(sqrt(1/sum(sr.weight)))
sr.combinedVarLogOR<-(1/sum(sr.weight))
sr.combinedChisq<-(((sr.combinedLogOR-0)^2)/sr.combinedVarLogOR)
sr.combinedValue<-pchisq(sr.combinedChisq, df=1)
sr.combinedPvalue<-(1-sr.combinedValue)

#get qnorm values
alpha<-(1-((1-catmapobject$ci)/2))
quantilenorm<-qnorm(alpha, 0, 1)

sr.lbci<-exp(sr.combinedLogOR-(quantilenorm*sr.combinedSeLogOR))
sr.ubci<-exp(sr.combinedLogOR+(quantilenorm*sr.combinedSeLogOR))
sr.combinedCI<-c(sr.lbci, sr.ubci)
sr.SeLogOR<-sqrt(sr.combinedVarLogOR)
sr.lbci<-exp(sr.logOR-(quantilenorm*sr.SeLogOR))
sr.ubci<-exp(sr.logOR+(quantilenorm*sr.SeLogOR))

#calculate heterogeneity
sr.df<-(nrow(catmapobject$a1)-2)
sr.chisqHet<-(sum(sr.weight*(((sr.logOR-sr.combinedLogOR)^2))))
sr.combinedHetValue<-pchisq(sr.chisqHet, df=sr.df)
sr.heterogeneityPvalue<-(1-sr.combinedHetValue)

#DerSimonian and Laird random-effects sensitivity analysis

sr.tau2<-((sr.chisqHet-sr.df)/(sum(sr.weight)-(sum(sr.weight^2)/(sum(sr.weight)))))

if (sr.tau2 <=0){
sr.tau2<-0
}
srweight.dsl<-(1/(sr.comvarlogOR+sr.tau2))
srlogOR.dsl<-((sum(srweight.dsl*sr.logOR))/(sum(srweight.dsl)))
srOR.dsl<-exp(srlogOR.dsl)
srseLogOR.dsl<-(1/(sqrt(sum(srweight.dsl))))
srvarLogOR.dsl<-(1/sum(srweight.dsl))
srlbci.dsl<-exp(srlogOR.dsl-(quantilenorm*srseLogOR.dsl))
srubci.dsl<-exp(srlogOR.dsl+(quantilenorm*srseLogOR.dsl))
srci.dsl<-c(srlbci.dsl, srubci.dsl)
srchisq.dsl<-(((srlogOR.dsl-0)^2)/srvarLogOR.dsl)
srvalue.dsl<-pchisq(srchisq.dsl, df=1)
srpvalue.dsl<-(1-srvalue.dsl)

srstudy.removed<-paste("Study Removed =", catmapobject$studyname[r], sep=" ")
srtable.header<-c("Q Statistic (Heterogeneity) Chi-Square", "Q Statistic (Heterogeneity) p-value", "DerSimonian & Laird Random-Effects OR", "DerSimonian & Laird Random-Effects Lower Bound CI", "DerSimonian & Laird Random-Effects Chi-Square", "DerSimonian & Laird Random-Effects p-value")
srtable.fill<-c(sr.chisqHet, sr.heterogeneityPvalue, srOR.dsl, srci.dsl, srpvalue.dsl)
sr.results<-rbind(srtable.header, round(srtable.fill, digits=5))
srvalues<-c(srOR.dsl, srci.dsl)
srplot[r,]<-srvalues
cat("Random Effects Sensitivity Analysis\n")
cat(srstudy.removed, sr.results, sep="\n")
cat("\n")

#print results to file dataset.random.effects.sensitivity.txt

if(catmapobject$dataset!=catmapdata){
sink(paste(catmapobject$dataset, "random.effects.sensitivity.txt", sep="."), append=TRUE)
cat("Random Effects Sensitivity Analysis\n")
cat(srstudy.removed, sr.results, sep="\n")
cat("\n")
sink()
}

if(catmapobject$dataset==catmapdata){
sink(paste("catmapdata", "random.effects.sensitivity.txt", sep="."), append=TRUE)
cat("Random Effects Sensitivity Analysis\n")
cat(srstudy.removed, sr.results, sep="\n")
cat("\n")
sink()
}
}
}

ci100<- catmapobject$ci*100

#create fixed-effects sensitivity plots
if (fe.sense==FALSE & fe.senseplot==TRUE){
cat("NOTICE: No fixed-efffects sensitivity analyses have been performed. \n Please set fe.sense==TRUE and try again.\n")
}

if (fe.senseplot==TRUE){
if(catmapobject$dataset!=catmapdata){
pdf(file=(paste(catmapobject$dataset, "fixed.effects.sensitivity.plot.pdf", sep=".")))
}
if(catmapobject$dataset==catmapdata){
pdf(file=(paste("catmapdata", "fixed.effects.sensitivity.plot.pdf", sep=".")))
}
sfpstudy<-c(1:nrow(catmapobject$a1))
sfplx<-max((min(sfplot[,2])-0.25),0)
sfphx<-max(sfplot[,3])+0.25
sfply<-min(sfpstudy)-0.5
sfphy<-max(sfpstudy)+0.5
mar1<-c(5.1, 7.1, 4.1, 2.1)
las1<-1
par("las"=las1)
par("mar"=mar1)
sfpdummy<-c(rep(NA, length(sfpstudy)))
sfpdummy[1]<- sfphx
sfpdummy[2]<- sfplx
sfpydummy<-c(rep(NA, length(sfpstudy)))
sfpydummy[1]<- sfphy
sfpydummy[2]<- sfply
xtitle<-paste("OR(", ci100, "% CI)")
plot(sfpdummy, sfpydummy, type="n", log="x", ylab="", ylim=rev(c(sfply, sfphy)), yaxt="n", xlab=xtitle, main="Sensitivity Analysis: \n Inverse Variance (Fixed-Effects) ORs")
abline(v=1.0)
for(z in 1:nrow(catmapobject$a1)){
points(sfplot[z,1],sfpstudy[z], pch=22, bg="black", col="black")
segments(sfplot[z,2],sfpstudy[z],sfplot[z,3],sfpstudy[z], bg="black", col="black")
sfpstudyname<-paste("Study Removed:", catmapobject$studyname[z], sep="\n")
mtext(paste(sfpstudyname),side=2, at=sfpstudy[z], cex=0.80)
}
#dev.off()
}

ci100<- catmapobject$ci*100

#create random-effects sensitivity plots
if (re.sense==FALSE & re.senseplot==TRUE){
cat("NOTICE: No random-efffects sensitivity analyses have been performed. \n Please set re.sense==TRUE and try again.\n")
}
if(fe.senseplot==TRUE){
dev.off()
}
if (re.senseplot==TRUE & catmapobject$tau2 > 0){
if(catmapobject$dataset!=catmapdata){
pdf(file=(paste(catmapobject$dataset, "random.effects.sensitivity.plot.pdf", sep=".")))
}
if(catmapobject$dataset==catmapdata){
pdf(file=(paste("catmapdata", "random.effects.sensitivity.plot.pdf", sep=".")))
}
srpstudy<-c(1:nrow(catmapobject$a1))
srplx<-max((min(srplot[,2])-0.25),0)
srphx<-max(srplot[,3])+0.25
srply<-min(srpstudy)-0.5
srphy<-max(srpstudy)+0.5
mar1<-c(5.1, 7.1, 4.1, 2.1)
las1<-1
par("las"=las1)
par("mar"=mar1)
srpdummy<-c(rep(NA, length(srpstudy)))
srpdummy[1]<- srphx
srpdummy[2]<- srplx
srpydummy<-c(rep(NA, length(srpstudy)))
srpydummy[1]<- srphy
srpydummy[2]<- srply
xtitle<-paste("OR(", ci100, "% CI)")
plot(srpdummy, srpydummy, type="n", log="x", ylab="", ylim=rev(c(srply, srphy)), yaxt="n", xlab=xtitle, main="Sensitivity Analysis: \n DerSimonian & Laird (Random-Effects) ORs")
abline(v=1.0)
for(y in 1:nrow(catmapobject$a1)){
points(srplot[y,1],srpstudy[y], pch=22, bg="black", col="black")
segments(srplot[y,2],srpstudy[y],srplot[y,3],srpstudy[y], bg="black", col="black")
srpstudyname<-paste("Study Removed:", catmapobject$studyname[y], sep="\n")
mtext(paste(srpstudyname),side=2, at=srpstudy[y], cex=0.80)
}
graphics.off()
}
}

catmap.cumulative<-function(catmapobject, fe.cumulative, re.cumulative, fe.cumplot, re.cumplot){

#fixed-effect cumulative meta analyses
ci100<- catmapobject$ci*100

if (fe.cumulative==TRUE){

cfplot<-matrix(0, nrow(catmapobject$a1),3)
for(c in 1:nrow(catmapobject$a1)){
cf.weight<- catmapobject$weight[1:c]
cf.logOR<- catmapobject$logOR[1:c]

cf.combinedLogOR<-((sum(cf.weight*cf.logOR))/sum(cf.weight))
cf.combinedOR<-exp(cf.combinedLogOR)
cf.combinedSeLogOR<-(sqrt(1/sum(cf.weight)))
cf.combinedVarLogOR<-(1/sum(cf.weight))
cf.combinedChisq<-(((cf.combinedLogOR-0)^2)/cf.combinedVarLogOR)
cf.combinedValue<-pchisq(cf.combinedChisq, df=1)
cf.combinedPvalue<-(1-cf.combinedValue)

#get qnorm values
alpha<-(1-((1-catmapobject$ci)/2))
quantilenorm<-qnorm(alpha, 0, 1)

cf.lbci<-exp(cf.combinedLogOR-(quantilenorm*cf.combinedSeLogOR))
cf.ubci<-exp(cf.combinedLogOR+(quantilenorm*cf.combinedSeLogOR))
cf.combinedCI<-c(cf.lbci, cf.ubci)
cf.SeLogOR<-sqrt(cf.combinedVarLogOR)
cf.lbci<-exp(cf.logOR-(quantilenorm*cf.SeLogOR))
cf.ubci<-exp(cf.logOR+(quantilenorm*cf.SeLogOR))

#calculate heterogeneity
cf.df<-(c-1)
cf.chisqHet<-(sum(cf.weight*(((cf.logOR-cf.combinedLogOR)^2))))
cf.combinedHetValue<-pchisq(cf.chisqHet, df=cf.df)
cf.heterogeneityPvalue<-(1-cf.combinedHetValue)

study.added<-paste("Study Added =", catmapobject$studyname[c], sep=" ")
cftable.header<-c("Inverse Variance Fixed-Effects OR", "Inverse Variance Fixed-Effects Lower Bound CI", "Inverse Variance Fixed-Effects Upper Bound CI", "Inverse Variance Fixed-Effects Chi-Square", "Inverse Variance Fixed-Effects p-value", "Q Statistic (Heterogeneity) Chi-Square", "Q Statistic (Heterogeneity) p-value")
cftable.fill<-c(cf.combinedOR, cf.combinedCI, cf.combinedChisq, cf.combinedPvalue, cf.chisqHet, cf.heterogeneityPvalue)
cf.results<-rbind(cftable.header, round(cftable.fill, digits=5))
cfvalues<-c(cf.combinedOR, cf.combinedCI)
cfplot[c,]<-cfvalues
cat("Fixed Effects Cumulative Meta-Analysis\n")
cat(study.added, cf.results, sep="\n")
cat("\n")

#print results to file dataset.fixed.effects.cumulative.txt

if(catmapobject$dataset!=catmapdata){
sink(paste(catmapobject$dataset, "fixed.effects.cumulative.txt", sep="."), append=TRUE)
cat("Fixed Effects Cumulative Meta-Analysis\n")
cat(study.added, cf.results, sep="\n")
cat("\n")
sink()
}

if(catmapobject$dataset==catmapdata){
sink(paste("catmapdata", "fixed.effects.cumulative.txt", sep="."), append=TRUE)
cat("Fixed Effects Cumulative Meta-Analysis\n")
cat(study.added, cf.results, sep="\n")
cat("\n")
sink()
}
}
}

#random-effect cumulative meta analyses

if(re.cumulative==TRUE & catmapobject$tau2 <= 0){
cat("NOTICE: tau2 is less than or equal to 0;\n no random effects estimates will be calculated.\n")
}

if (re.cumulative==TRUE & catmapobject$tau2 > 0){
crplot<-matrix(0,(nrow(catmapobject$a1)-1),3)
for(u in 2:nrow(catmapobject$a1)){
v<-u-1
cr.weight<- catmapobject$weight[1:u]
cr.logOR<- catmapobject$logOR[1:u]
cr.comvarlogOR<- catmapobject$comvarlogOR[1:u]

cr.combinedLogOR<-((sum(cr.weight*cr.logOR))/sum(cr.weight))
cr.combinedOR<-exp(cr.combinedLogOR)
cr.combinedSeLogOR<-(sqrt(1/sum(cr.weight)))
cr.combinedVarLogOR<-(1/sum(cr.weight))
cr.combinedChisq<-(((cr.combinedLogOR-0)^2)/cr.combinedVarLogOR)
cr.combinedValue<-pchisq(cr.combinedChisq, df=1)
cr.combinedPvalue<-(1-cr.combinedValue)

#get qnorm values
alpha<-(1-((1- catmapobject$ci)/2))
quantilenorm<-qnorm(alpha, 0, 1)

cr.lbci<-exp(cr.combinedLogOR-(quantilenorm*cr.combinedSeLogOR))
cr.ubci<-exp(cr.combinedLogOR+(quantilenorm*cr.combinedSeLogOR))
cr.combinedCI<-c(cr.lbci, cr.ubci)
cr.SeLogOR<-sqrt(cr.combinedVarLogOR)

#calculate heterogeneity
cr.df<-(u-1)
cr.chisqHet<-(sum(cr.weight*(((cr.logOR-cr.combinedLogOR)^2))))
cr.combinedHetValue<-pchisq(cr.chisqHet, df=cr.df)
cr.heterogeneityPvalue<-(1-cr.combinedHetValue)

#DerSimonian and Laird random-effects cumulative analysis

cr.tau2c<-(cr.chisqHet-cr.df)/(sum(cr.weight)-(sum(cr.weight^2)/(sum(cr.weight))))

cr.tau2<-max(0,cr.tau2c)

if (cr.tau2 <=0){
cr.tau2<-0
}

crweight.dsl<-(1/(cr.comvarlogOR+cr.tau2))
crlogOR.dsl<-((sum(crweight.dsl*cr.logOR))/(sum(crweight.dsl)))
crOR.dsl<-exp(crlogOR.dsl)
crseLogOR.dsl<-(1/(sqrt(sum(crweight.dsl))))
crvarLogOR.dsl<-(1/sum(crweight.dsl))
crlbci.dsl<-exp(crlogOR.dsl-(quantilenorm*crseLogOR.dsl))
crubci.dsl<-exp(crlogOR.dsl+(quantilenorm*crseLogOR.dsl))
crci.dsl<-c(crlbci.dsl, crubci.dsl)
crchisq.dsl<-(((crlogOR.dsl-0)^2)/crvarLogOR.dsl)
crvalue.dsl<-pchisq(crchisq.dsl, df=1)
crpvalue.dsl<-(1-crvalue.dsl)

crstudy.added<-paste("Study Added =", catmapobject$studyname[u], sep=" ")
crtable.header<-c("Q Statistic (Heterogeneity) Chi-Square", "Q Statistic (Heterogeneity) p-value", "DerSimonian & Laird Random-Effects OR", "DerSimonian & Laird Random-Effects Lower Bound CI", "DerSimonian & Laird Random-Effects Chi-Square", "DerSimonian & Laird Random-Effects p-value")
crtable.fill<-c(cr.chisqHet, cr.heterogeneityPvalue, crOR.dsl, crci.dsl, crpvalue.dsl)
cr.results<-rbind(crtable.header, round(crtable.fill, digits=5))
crvalues<-c(crOR.dsl, crci.dsl)
crplot[v,]<-crvalues
cat("Random Effects Cumulative Meta-Analysis\n")
cat(crstudy.added, cr.results, sep="\n")
cat("\n")

#print results to file dataset.random.effects.cumulative.txt

if(catmapobject$dataset!=catmapdata){
sink(paste(catmapobject$dataset, "random.effects.cumulative.txt", sep="."), append=TRUE)
cat("Random Effects Cumulative Meta-Analysis\n")
cat(crstudy.added, cr.results, sep="\n")
cat("\n")
sink()
}

if(catmapobject$dataset==catmapdata){
sink(paste("catmapdata", "random.effects.cumulative.txt", sep="."), append=TRUE)
cat("Random Effects Cumulative Meta-Analysis\n")
cat(crstudy.added, cr.results, sep="\n")
cat("\n")
sink()
}
}
}

#create fixed-effects cumulative plots

if(fe.cumulative==FALSE & fe.cumplot==TRUE){
cat("NOTICE: No fixed-efffects cumulative analyses have been performed.\n Please set fe.cumulative==TRUE and try again.\n")
}

if (fe.cumplot==TRUE){
if(catmapobject$dataset!=catmapdata){
pdf(file=(paste(catmapobject$dataset, "fixed.effects.cumulative.plot.pdf", sep=".")))
}
if(catmapobject$dataset==catmapdata){
pdf(file=(paste("catmapdata", "fixed.effects.cumulative.plot.pdf", sep=".")))
}
cfpstudy<-c(1:nrow(catmapobject$a1))
cfplx<-max((min(cfplot[,2])-0.25),0)
cfphx<-max(cfplot[,3])+0.25
cfply<-1-0.5
cfphy<-max(cfpstudy)+0.5
mar1<-c(5.1, 7.1, 4.1, 2.1)
las1<-1
par("las"=las1)
par("mar"=mar1)
cfpdummy<-c(rep(NA, length(cfpstudy)))
cfpdummy[1]<- cfphx
cfpdummy[2]<- cfplx
cfpydummy<-c(rep(NA, length(cfpstudy)))
cfpydummy[1]<- cfphy
cfpydummy[2]<- cfply
xtitle<-paste("OR(", ci100, "% CI)")
plot(cfpdummy, cfpydummy, type="n", log="x", ylab="", ylim=rev(c(cfply, cfphy)), yaxt="n", xlab=xtitle, main="Cumulative Meta-Analysis: \n Inverse Variance (Fixed-Effects) ORs")
abline(v=1.0)
for(x in 1:nrow(catmapobject$a1)){
points(cfplot[x,1],cfpstudy[x], pch=22, bg="black", col="black")
segments(cfplot[x,2],cfpstudy[x],cfplot[x,3],cfpstudy[x], bg="black", col="black")
cfpstudyname<-paste("Study Added:", catmapobject$studyname[x], sep="\n")
mtext(paste(cfpstudyname),side=2, at=cfpstudy[x], cex=0.80)
}
dev.off()
}

#create random-effects cumulative plots
if(re.cumulative==FALSE & re.cumplot==TRUE){
cat("NOTICE: No random-efffects cumulative analyses have been performed. \n Please set re.cumulative==TRUE and try again.\n")
}
if (re.cumplot==TRUE & catmapobject$tau2 > 0){
if(catmapobject$dataset!=catmapdata){
pdf(file=(paste(catmapobject$dataset, "random.effects.cumulative.plot.pdf", sep=".")))
}
if(catmapobject$dataset==catmapdata){
pdf(file=(paste("catmapdata", "random.effects.cumulative.plot.pdf", sep=".")))
}
crpstudy<-c(1:(nrow(catmapobject$a1)-1))
crplx<-max((min(crplot[,2])-0.25),0)
crphx<-max(crplot[,3])+0.25
crply<-min(crpstudy)-0.5
crphy<-max(crpstudy)+0.5
mar1<-c(5.1, 7.1, 4.1, 2.1)
las1<-1
par("las"=las1)
par("mar"=mar1)
crpdummy<-c(rep(NA, length(crpstudy)))
crpdummy[1]<- crphx
crpdummy[2]<- crplx
crpydummy<-c(rep(NA, length(crpstudy)))
crpydummy[1]<- crphy
crpydummy[2]<- crply
xtitle<-paste("OR (", ci100, "% CI)")
plot(crpdummy, crpydummy, type="n", log="x", ylab="", ylim=rev(c(crply, crphy)), yaxt="n", xlab=xtitle, main="Cumulative Meta-Analysis: \n DerSimonian & Laird (Random-Effects) ORs")
abline(v=1.0)
for(w in 1:(nrow(catmapobject$a1)-1)){
q<-w+1
points(crplot[w,1],crpstudy[w], pch=22, bg="black", col="black")
segments(crplot[w,2],crpstudy[w],crplot[w,3],crpstudy[w], bg="black", col="black")
crpstudyname<-paste("Study Added:", catmapobject$studyname[w], sep="\n")
mtext(paste(crpstudyname),side=2, at=crpstudy[w], cex=0.80)
}
graphics.off()
}
}

# create fixed-effect forest plots

catmap.forest<-function(catmapobject, fe.forest, re.forest){
ci100<- catmapobject$ci*100
if (fe.forest==TRUE){
if(catmapobject$dataset!=catmapdata){
pdf(file=(paste(catmapobject$dataset, "fixed.effects.plot.pdf", sep=".")))
}
if(catmapobject$dataset==catmapdata){
pdf(file=(paste("catmapdata", "fixed.effects.plot.pdf", sep=".")))
}
prop.weight<-(catmapobject$weight/(sum(catmapobject$weight))*10)
OR<-c(catmapobject$combinedOR, exp(catmapobject$logOR))
study<-1:(nrow(catmapobject$a1)+1)
lx<-max((min(catmapobject$lbci.fe)-0.25),0)
hx<-max(catmapobject$ubci.fe)+0.25
ly<-(1-0.5)
hy<-(length(catmapobject$study)+ 1.5)
mar1<-c(5.1, 7.1, 4.1, 2.1)
las1<-1
par("las"=las1)
par("mar"=mar1)
dummy<-c(rep(NA, length(catmapobject$study)))
dummy[1]<-hx
dummy[2]<-lx
ydummy<-c(rep(NA, length(catmapobject$study)))
ydummy[1]<-hy
ydummy[2]<-ly
xtitle<-paste("OR(", ci100, "% CI)")
plot(dummy, ydummy, type="n", log="x", ylab="", yaxt="n", xlab=xtitle, main="Inverse Variance (Fixed-Effects) ORs")
points(OR[1], 1, pch=22, cex=1.0, bg="red", col="red")
segments(catmapobject$lbci, 1, catmapobject$ubci, 1, col="red")
mtext("Pooled OR",side=2,at=1, cex=0.80)
abline(v=1.0)
limit<-nrow(catmapobject$a1)+1
counts<-c(1:limit)
for(i in 2:(nrow(catmapobject$a1)+1)){
j<-i-1
points(OR[i], counts[i], pch=22, cex=prop.weight[j], bg="black", col="black")
segments(catmapobject$lbci.fe[j], counts[i], catmapobject$ubci.fe[j], counts[i], bg="black", col="black")
mtext(paste(catmapobject$studyname[j]),side=2, at=counts[i], cex=0.80)
}
dev.off()
}

#create random-effects forest plots 

if(re.forest==TRUE & catmapobject$tau2 <= 0){
cat("NOTICE: tau2 is less than or equal to 0;\n no random effects estimates will be calculated.\n")
}

if (re.forest==TRUE & catmapobject$tau2 > 0){
if(catmapobject$dataset!=catmapdata){
pdf(file=(paste(catmapobject$dataset, "random.effects.plot.pdf", sep=".")))
}
if(catmapobject$dataset==catmapdata){
pdf(file=(paste("catmapdata", "random.effects.plot.pdf", sep=".")))
}
dslprop.weight<-(catmapobject$weight.dsl/(sum(catmapobject$weight.dsl))*10)
dslOR<-c(catmapobject$OR.dsl, exp(catmapobject$logOR))
OR<-c(catmapobject$combinedOR, exp(catmapobject$logOR))
rstudy<-1:(nrow(catmapobject$a1)+1)
rlx<-max((min(catmapobject$lbci.fe)-0.25),0)
rhx<-max(catmapobject$ubci.fe)+0.25
rly<-(1-0.5)
rhy<-(length(catmapobject$study)+1.5)
mar1<-c(5.1, 7.1, 4.1, 2.1)
las1<-1
par("las"=las1)
par("mar"=mar1)
rdummy<-c(rep(NA, length(rstudy)))
rdummy[1]<-rhx
rdummy[2]<-rlx
rydummy<-c(rep(NA, length(rstudy)))
rydummy[1]<-rhy
rydummy[2]<-rly
xtitle<-paste("OR(", ci100, "% CI)")
plot(rdummy, rydummy, type="n", log="x", ylab="", yaxt="n", xlab=xtitle, main="DerSimonian & Laird (Random-Effects) ORs")
points(dslOR[1], 1, pch=22, cex=1.0, bg="red", col="red")
segments(catmapobject$lbci.dsl, 1, catmapobject$ubci.dsl, 1, col="red")
mtext("Pooled OR",side=2,at=1, cex=0.80)
abline(v=1.0)
limit<-nrow(catmapobject$a1)+1
counts<-c(1:limit)
for(e in 2:(nrow(catmapobject$a1)+1)){
f<-e-1
points(OR[e],counts[e], pch=22, cex=dslprop.weight[f], bg="black", col="black")
segments(catmapobject$lbci.fe[f],counts[e], catmapobject$ubci.fe[f],counts[e], bg="black", col="black")
mtext(paste(catmapobject$studyname[f]),side=2, at=counts[e], cex=0.80)
}
graphics.off()
}
}

#funnel plot

catmap.funnel<-function(catmapobject, funnel){
if(funnel==TRUE){
if(catmapobject$dataset!=catmapdata){
pdf(file=(paste(catmapobject$dataset, "funnel.plot.pdf", sep=".")))
}
if(catmapobject$dataset==catmapdata){
pdf(file=(paste("catmapdata", "funnel.plot.pdf", sep=".")))
}
OddsRatio<-exp(catmapobject$logOR)
fun.se<-sqrt(catmapobject$comvarlogOR)
maxse<-max(fun.se)
minse<-min(fun.se)
plot(OddsRatio, fun.se, ylim=rev(c(minse,maxse)),log="x", ylab="Standard Error Log OR")
abline(v=1)
abline(v=catmapobject$combinedOR, lty=2)
graphics.off()
}
}
