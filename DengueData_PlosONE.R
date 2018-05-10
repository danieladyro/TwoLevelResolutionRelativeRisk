rm(list=ls())

################################################
## To download and install the R-INLA package ##
##    (http://www.r-inla.org/download)        ##
################################################
# install.packages("INLA", repos=c(getOption("repos"),
#                  INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(INLA)
library(ggplot2)
library(rgdal)
library(broom)
library(gridExtra)
library(car)

## n=94 census sector
## m=17 communes
## t=91 epidemiological periods
n <- 94
m <- 17
t <- 91
t.from <- 2009
t.to <- 2015


########################################################
## Spatial neighborhood matrix for the CENSUS SECTORS ##
########################################################
g <- inla.read.graph("sectoresGraph.dat")
R.xi <- matrix(0, g$n, g$n)
for (i in 1:g$n){
	R.xi[i,i]=g$nnbs[[i]]
	R.xi[i,g$nbs[[i]]]=-1
}

## Spatial structure matrix for a LCAR prior ##
R.Leroux.FLA <- diag(n)-R.xi


##################################################
## Spatial neighborhood matrix for the COMMUNES ##
##################################################
g <- inla.read.graph("comunasGraph.dat")
R.psi <- matrix(0, g$n, g$n)
for (i in 1:g$n){
	R.psi[i,i]=g$nnbs[[i]]
	R.psi[i,g$nbs[[i]]]=-1
}

## Spatial structure matrix for a LCAR prior ##
R.Leroux.SLA <- diag(m)-R.psi


###############################################
## Temporal structure matrix for a RW1 prior ##
###############################################
D1 <- diff(diag(t),differences=1)
R.gammaRW1 <- t(D1)%*%D1


########################################################
## Cases of dengue disease from Bucaramanga, Colombia ##
########################################################
load("DengueData.RData")

## "DengueData_INLA" dataframe variables:
## O = Number of observed cases (per census sector and epidemiological period)
## E = Number of expected cases (per census sector and epidemiological period)
## ID.FLA = IDs for the census sectors (first-level area)
## ID.SLA = IDs for the communes (second-level area)
## ID.year = IDs for the epidemiological periods
## ID.FLA.year = IDs for the space-time interaction
## sector = Original IDs (shapefile) for the census sectors


##################################################################################################
## Figure 1: Descriptive analysis of dengue disease cases in the city of Bucaramanga, Colombia. ##
##################################################################################################

## (A) Cases by epidemiological period.
temporaltotal <- aggregate(DengueData_INLA$O,list(DengueData_INLA$ID.year),sum)
colnames(temporaltotal) <- c("Time","cases")

theme_longitudinal <-	theme_bw() + theme(
  panel.border = element_rect(colour="black", size=0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=18, margin=margin(15,5,10,5,"pt")),  
  axis.text.x  = element_text(size=14), 
  axis.title.y = element_text(size=18, margin=margin(10,15,10,5,"pt")),
  axis.text.y  = element_text(size=14),
  plot.title   = element_text(hjust=-0.2, size=24),
  legend.text  = element_text(size=14,lineheight=1.0),
  legend.title = element_text(size=14),
  plot.caption = element_text(hjust=0.5, size=14),
  legend.key.height=unit(1, "cm"),
  strip.text = element_text(face="bold", size=14),
  strip.background = element_rect(fill="white", colour="white")) 

temporalall <- ggplot(temporaltotal, aes(Time,cases)) + 
  geom_line(colour="blue",size=1) +
  scale_x_continuous(breaks=c(1,14,27,40,53,66,79),
                     labels=c(2009,2010,2011,2012,2013,2014,2015)) +
  ylab("Cases") +
  xlab("Epidemiological year")+
  labs(title="(A)") +
  theme_longitudinal


## (B) Annual average cases per 100.000 inhabitants and age-groups.
datacuminc <- data.frame(cuminc=cuminc$X2015[1:14]/7, age=seq(1,14,1))

theme_bar <-	theme_bw() + theme(
  panel.border = element_rect(colour="black", size=0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=18, margin=margin(15,5,10,5,"pt")),  
  axis.text.x  = element_text(size=14,angle=45,vjust=0.6), 
  axis.title.y = element_text(size=18, margin=margin(10,15,10,5,"pt")),
  axis.text.y  = element_text(size=14),
  plot.title   = element_text(hjust=-0.25, size=24),
  legend.text  = element_text(size=14,lineheight=1.0),
  legend.title = element_text(size=14),
  plot.caption = element_text(hjust=0.5, size=14),
  legend.key.height=unit(1, "cm"),
  strip.text = element_text(face="bold", size=14),
  strip.background = element_rect(fill="white", colour="white")) 

cummulative <- ggplot(data=datacuminc,aes(x=age,y=cuminc)) + 
  geom_bar(colour="white",fill="blue",stat="identity")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),
                     labels=c("0-4 ","5-9 ","10-14 ","15-19 ","20-25 ","24-29 ","30-34 ",
                              "35-39 ","40-44 ","45-49 ","50-54 ","55-59 ","60-64 ","65+ ")) +
  ylab("Annual average cases per 100.000 inhabitants") +
  xlab("Age-group")+
  labs(title="(B)") +
  theme_bar

postscript("Fig1.eps", width=12, height=9)
grid.arrange(temporalall,cummulative,nrow=1)
dev.off()


################################################################################################################
## Figure 2: Cumulative standardized incidence rates (SIRs) of dengue disease by communes and census sectors. ##
################################################################################################################
Carto.sector.np <- readOGR("sector293NP/.", "sector293NP")
map.sector.np <- tidy(Carto.sector.np, region="ID")

Carto.comuna.np <- readOGR("comunaNP/.", "comunaNP")
Carto.comuna.np@data$ID.comuna <- seq(1,17,1)
map.comuna.np <- tidy(Carto.comuna.np, region ="ID.comuna")


## (A) SIR of dengue cases by census sector.
casesBySector <- data.frame(sector=unique(DengueData_INLA$sector),
                            casos=aggregate(DengueData_INLA$O,list(DengueData_INLA$ID.FLA),sum)$x,
                            expcasos=aggregate(DengueData_INLA$E,list(DengueData_INLA$ID.FLA),sum)$x)

sectorDengue <- merge(Carto.sector.np@data, casesBySector, by.x="ID", by.y="sector")

holes <- geom_map(#inherit.aes = FALSE, #fill the holes with white
  aes(map_id=id),
  data = map.sector.np[which(map.sector.np$hole==TRUE),],
  map  = map.sector.np[which(map.sector.np$hole==TRUE),],
  fill = "#FFFFFF",
  alpha = 1)

mytheme <- 	theme_bw() + theme(
  panel.border = element_rect(colour="black", size=0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x =element_text(size=18),  
  axis.text.x  = element_text(size=18), 
  axis.title.y =element_text(size=18),
  axis.text.y  = element_text(size=18),
  plot.title   = element_text(hjust=-0.25, size=24),
  legend.text  =element_text(size=18,lineheight=0.8),
  legend.title =element_text(size=18),
  legend.key.height=unit(1, "cm"),
  plot.caption = element_text(hjust=0.5, size=22),
  strip.text = element_text(face="bold", size=18),
  strip.background = element_rect(fill="white", colour="white")) 

mapSIRSector <- ggplot() + 
  geom_map(data = sectorDengue, aes(map_id=ID, fill=casos/expcasos),
           map=map.sector.np,col="black",size=0.05) +
  expand_limits(x = map.sector.np$long, y = map.sector.np$lat) + 
  scale_fill_continuous(low='aquamarine',high='red') +
  scale_x_continuous(breaks=c(-73.15,-73.11),
                     labels=c(-73.15, -73.11)) +
  labs(fill="SIR",title="(A)") + 
  holes +
  ylab("Longitude") + 
  xlab("Latitude")	+
  mytheme


## (B) SIR of dengue cases by commune.
Data2Com <- merge(DengueData_INLA, Carto.sector.np@data, by.x="sector", by.y="ID")

casesBycomuna <- aggregate(list(Data2Com$O,Data2Com$E), by=list(ID=Data2Com$comuna), sum)
casesBycomuna$ID.comuna <- seq(1,17,1)
names(casesBycomuna) <- c("comuna","casos","expcasos","ID.comuna")

holes.com <- geom_map(#inherit.aes = FALSE, #fill the holes with white
  aes(map_id=id),
  data = map.comuna.np[which(map.comuna.np$hole==TRUE),],
  map=   map.comuna.np[which(map.comuna.np$hole==TRUE),],
  fill = "#FFFFFF",
  alpha = 1)

mapSIRComuna <- ggplot() + 
  geom_map(data=casesBycomuna, aes(map_id=ID.comuna, fill=casos/expcasos),
           map=map.comuna.np,size=0.05,col="black") +
  expand_limits(x = map.comuna.np$long, y = map.comuna.np$lat) + 
  scale_fill_continuous(low='aquamarine',high='red') +
  scale_x_continuous(breaks=c(-73.15,-73.11),
                     labels=c(-73.15, -73.11)) +
  # labs(fill="Cases",caption="SIR of dengue by commune") +
  labs(fill="SIR",title="(B)") +
  holes.com+
  ylab("Longitude") + 
  xlab("Latitude")	+
  mytheme


postscript("Fig2.eps", width=12, height=9)
grid.arrange(mapSIRSector,mapSIRComuna,nrow=1)
dev.off()


#########################################
## Define the hyperprior distributions ##
#########################################
sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

lunif = "expression:
    a = 1;
    b = 1;
    beta = exp(theta)/(1+exp(theta));
    logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
    log_jacobian = log(beta*(1-beta));
    return(logdens+log_jacobian)"


########################################################################
######		               Two-level model A		                    ######
######	log(r_it) = alpha + csi_i + psi_j(i) + gamma_t + delta_it	######
########################################################################
R <- kronecker(R.xi,R.gammaRW1)
r.def <- n+t-1
A1 <- kronecker(diag(n),matrix(1,1,t))
A2 <- kronecker(matrix(1,1,n),diag(t))
delta.constr <- rbind(A1,A2)

formula.TLA.IV <- O ~ f(ID.FLA, model="generic1", Cmatrix = R.Leroux.FLA, constr=TRUE,
                        hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
                      f(ID.SLA, model="generic1", Cmatrix = R.Leroux.SLA, constr=TRUE,
                        hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
                      f(ID.year, model="rw1", hyper=list(prec=list(prior=sdunif))) +
                      f(ID.FLA.year, model="generic0", Cmatrix=R, rankdef=r.def,
                        constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                        extraconstr=list(A=delta.constr, e=rep(0,n+t)))

Model.TLA.IV <- inla(formula.TLA.IV, family="poisson", data=DengueData_INLA, E=E,
                     control.predictor=list(compute=TRUE, cdf=c(log(1))),
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     control.inla=list(strategy="simplified.laplace"))

## NOTE: To speed up computations -> strategy="gaussian" ##
## NOTE: Due to file size restrictions in the data repository,
##       the INLA object Model.TLA.IV has been edited, to keep only
##       the objects employed to make the paper's plots and tables

##############################################################################################
###############################################################################################
###############################################################################################

load("Dengue_INLA_TLModelA_TypeIV_edited.RData")

## Results of Table 2 ##
data.frame(D=Model.TLA.IV$dic$mean.deviance,
           pD=Model.TLA.IV$dic$p.eff,
           DIC=Model.TLA.IV$dic$dic,
           DICc=Model.TLA.IV$dic$mean.deviance+sum(Model.TLA.IV$dic$local.p.eff/(1-Model.TLA.IV$dic$local.p.eff)),
           WAIC=Model.TLA.IV$waic$waic,
           LS=-sum(log(Model.TLA.IV$cpo$cpo)))


## Results of Table 3 ##
data.frame(Model.TLA.IV$summary.hyperpar[,c(1:5)])


###########################################################################################################
## Figure 3: Posterior mean estimates of spatial random effects at both census sector and commune-level, ##
##           and posterior exceedance probability of being greater than one.                             ##
###########################################################################################################

## (A) Map of census sector level spatial incidence risk pattern exp(xi_i)
xi_exp <- unlist(lapply(Model.TLA.IV$marginals.random$ID.FLA, function(x) inla.emarginal(exp,x)))

datosspatpatt <- data.frame(region=unique(DengueData_INLA$ID.FLA), dat=xi_exp, sector=unique(DengueData_INLA$sector))
datosspatpatt.np <- merge(Carto.sector.np@data, datosspatpatt, by.x="ID", by.y="sector")

themeAppendixSp <-	theme_bw() + theme(
  panel.border = element_rect(colour="black", size=0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x =element_text(size=14),  
  axis.text.x  = element_text(size=14), 
  axis.title.y = element_text(size=14),
  #axis.ticks = element_blank(),
  axis.text.y  = element_text(size=14),
  plot.title   = element_text(hjust=-0.25, size=24),
  legend.text  =element_text(size=12,lineheight=1.0),
  legend.title =element_text(size=14),
  legend.key.height=unit(1, "cm"),
  plot.caption = element_text(hjust=0.5, size=14),
  strip.text = element_text(face="bold", size=14),
  strip.background = element_rect(fill="white", colour="white")) 

FLA.mean <- ggplot() +
  geom_map(data=datosspatpatt.np, map=map.sector.np, aes(fill=dat, 
                                                         map_id=ID), col="black",size=0.1) +
  expand_limits(x= map.sector.np$long, y= map.sector.np$lat ) + 
  holes + 
  themeAppendixSp +
  scale_fill_gradient2(low='blue4',mid='white',high='red4',midpoint=1.0) +
  scale_x_continuous(breaks=c(-73.15,-73.11),
                     labels=c(-73.15, -73.11)) +
  ylab("Longitude") + 
  xlab("Latitude")+
  labs(fill=NULL,title="(A)") 


## (B) Posterior probability distribution P(exp(xi_i)>1 | O)
prob.xi <- unlist(lapply(Model.TLA.IV$marginals.random$ID.FLA, function(x){1-inla.pmarginal(0, x)}))

datosprobspatpatt <- data.frame(region=unique(DengueData_INLA$ID.FLA), dat=prob.xi, sector=unique(DengueData_INLA$sector))
datosprobspatpatt.np <- merge(Carto.sector.np@data, datosprobspatpatt, by.x="ID", by.y="sector")

FLA.prob <- ggplot() +
  geom_map(data=datosprobspatpatt.np,  
           aes(map_id=ID,fill=dat), map=map.sector.np, col="black",size=0.1) +  
  expand_limits(x=map.sector.np$long,y=map.sector.np$lat) +
  holes+
  themeAppendixSp +
  scale_fill_gradient2(low='green4',mid='white',high='red4',midpoint=0.5,limits=c(0,1)) +
  scale_x_continuous(breaks=c(-73.15,-73.11),
                     labels =c(-73.15, -73.11)) +
  ylab("Longitude") + xlab("Latitude")+
  labs(fill=NULL,title="(B)") 


## (C) Map of commune level spatial incidence risk pattern exp(psi_j(i)).
psi_exp <- unlist(lapply(Model.TLA.IV$marginals.random$ID.SLA, function(x) inla.emarginal(exp,x)))

datosspcommune <- data.frame(id=factor(seq(1,17,1)), dat=psi_exp)
datosspcommune.np <- merge(Carto.comuna.np@data, datosspcommune, by.x="ID.comuna", by.y="id")

SLA.mean <- ggplot() +
  geom_map(data=datosspcommune.np, map=map.comuna.np, aes(fill=dat, map_id=ID.comuna), col="black",size=0.1) + 
  expand_limits(x=map.comuna.np$long,y=map.comuna.np$lat) +
  holes.com +
  themeAppendixSp +
  scale_fill_gradient2(low='blue4',mid='white',high='red4',midpoint=1.0) +
  scale_x_continuous(breaks=c(-73.15,-73.11),
                     labels=c(-73.15, -73.11)) +
  ylab("Longitude") + xlab("Latitude")+
  labs(fill=NULL,title="(C)") 


## (D) Posterior probability distribution P(exp(psi_j(i))>1 | O)
prob.psi <- unlist(lapply(Model.TLA.IV$marginals.random$ID.SLA, function(x){1-inla.pmarginal(0,x)}))

datosProbCom <- data.frame(id=factor(seq(1,17,1)), dat=prob.psi)
datosProbCom.np <- merge(Carto.comuna.np@data, datosProbCom, by.x="ID.comuna", by.y="id")

SLA.prob <- ggplot() +
  geom_map(data=datosProbCom.np, map=map.comuna.np, aes(fill=dat, map_id=ID.comuna), col="black",size=0.1) +
  expand_limits(x=map.comuna.np$long,y=map.comuna.np$lat) +
  holes.com +
  themeAppendixSp +
  scale_fill_gradient2(low='green4',mid='white',high='red4',midpoint=0.5,limits=c(0,1)) +
  scale_x_continuous(breaks=c(-73.15,-73.11), labels=c(-73.15, -73.11)) +
  ylab("Longitude") + 
  xlab("Latitude")+
  labs(fill=NULL, title="(D)") 


postscript("Fig3.eps", width=8)
grid.arrange(FLA.mean,FLA.prob,SLA.mean,SLA.prob,nrow=2)
dev.off()


##########################################################################################################
## Figure 4: Overall temporal trend of dengue disease incidence relative risk by epidemiological period ##
##########################################################################################################
temporal <- unlist(lapply(Model.TLA.IV$marginals.random$ID.year, function(x) inla.emarginal(exp,x)))
aux <- lapply(Model.TLA.IV$marginals.random$ID.year, function(x) inla.tmarginal(exp,x))
q1 <- unlist(lapply(aux, function(x) inla.qmarginal(0.025,x)))
q2 <- unlist(lapply(aux, function(x) inla.qmarginal(0.975,x)))

pl = data.frame(Time=1:91, menle=temporal)
pl$menlelb = q1
pl$menleub = q2

theme_longitudinal <-	theme_bw() + theme(
  panel.border = element_rect(colour="black", size=0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=26, margin=margin(15,5,10,5,"pt")),
  axis.text.x  = element_text(size=20), 
  axis.title.y = element_text(size=26),
  axis.text.y  = element_text(size=20),
  plot.title   = element_text(size=16),
  legend.text  = element_text(size=16,lineheight=1.0),
  legend.title = element_text(size=16),
  legend.key.height=unit(1, "cm"),
  plot.caption = element_text(hjust=0.5, size=18),
  strip.text = element_text(face="bold", size=10),
  strip.background = element_rect(fill="white", colour="white")) 

temporal.overall <- ggplot(pl, aes(Time)) + 
  geom_line(aes(y=menle), colour="blue") + 
  geom_ribbon(aes(ymin=menlelb, ymax=menleub), alpha=0.2,color="grey")+
  scale_x_continuous(breaks=c(1,14,27,40,53,66,79),
                     labels=c(2009,2010,2011,2012,2013,2014,2015)) +
  ylab(expression(exp(gamma[t]))) + xlab("Epidemiological year")+
  geom_hline(yintercept = 1,linetype="dashed") + 
  theme_longitudinal

cairo_ps("Fig4.eps",height=8,width=11)
temporal.overall
dev.off()


#####################################################################################
## Figure 5: Maps with the estimated posterior mean values of the relative risk of ##
## dengue disease by census sector for the epidemiological periods 1 to 8 of 2013. ##
#####################################################################################
mtlA.IVfitted     <- data.frame(DengueData_INLA, Model.TLA.IV$summary.fitted.values)
mtlA.IVfitted$prp <- 1-Model.TLA.IV$summary.fitted.values[,"1 cdf"]

ep2013 <- c(53,54,55,56,57,58,59,60) # Epidemiological periods of 2013
mtlA.IVslp13short <- mtlA.IVfitted[mtlA.IVfitted$ID.year %in% ep2013,]
mtlA.IVslp13short$ID.yearLabel <- mtlA.IVslp13short$ID.year
mtlA.IVslp13short$ID.yearLabel <- recode(mtlA.IVslp13short$ID.year,
                                         "'53'= c('EP 01, Dec 30 - Jan 26');
                                         '54'= c('EP 02, Jan 27 - Feb 23');
                                         '55'= c('EP 03, Feb 24 - Mar 23');
                                         '56'= c('EP 04, Mar 24 - Apr 20');
                                         '57'= c('EP 05, Apr 21 - May 18');
                                         '58'= c('EP 06, May 19 - Jun 15');
                                         '59'= c('EP 07, Jun 16 - Jul 14');
                                         '60'= c('EP 08, Jul 14 - Aug 10')")

themeAppendixShort <-	theme_bw() + theme(
  panel.border = element_rect(colour="black", size=0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=10),  
  axis.text.x  = element_text(size=10), 
  axis.title.y = element_text(size=10),
  axis.ticks   = element_blank(),
  axis.text.y  = element_text(size=10),
  plot.title   = element_text(size=14),
  legend.text  = element_text(size=12,lineheight=1.0),
  legend.title = element_text(size=18),
  legend.key.height=unit(1, "cm"),
  strip.text = element_text(face="bold", size=10),
  strip.background = element_rect(fill="white", colour="white")) 

holes <- geom_map(#inherit.aes = FALSE, #fill the holes with white
  aes(map_id=id),
  data = map.sector.np[which(map.sector.np$hole==TRUE),],
  map  = map.sector.np[which(map.sector.np$hole==TRUE),],
  fill = "#FFFFFF",
  alpha = 1)

risk2013.map <- ggplot() +
  geom_map(data=mtlA.IVslp13short, map=map.sector.np,
           aes(fill=mean, map_id=sector),color="black",size=0.1) +  
  expand_limits(x=map.sector.np$long,y=map.sector.np$lat) +
  holes +
  facet_wrap(~ID.yearLabel,ncol=4) +
  themeAppendixShort + 
  scale_fill_gradient2(low='green4',mid='white',high='red4',midpoint=1) +
  scale_x_continuous(breaks=c(-73.15,-73.11),
                     labels=c(-73.15, -73.11)) +
  ylab("Longitude") + xlab("Latitude") +
  labs(fill=NULL) + ggtitle("2013") 


postscript("Fig5.eps",height=8,width=11)
print(risk2013.map)
dev.off()


#####################################################################################
## Figure 6: Maps of the posterior probability distribution P(r_it > 1 | O) of     ##
## dengue disease by census sector for the epidemiological periods 1 to 8 of 2013. ##
#####################################################################################
prob2013.map <- ggplot() +
  geom_map(data=mtlA.IVslp13short, map=map.sector.np,
           aes(fill=prp, map_id=sector),color="black",size=0.1) +
  expand_limits(x=map.sector.np$long,y=map.sector.np$lat) +
  holes +
  facet_wrap(~ID.yearLabel,ncol=4) +
  themeAppendixShort + 
  scale_fill_gradient2(low='cyan',mid='white',high='red4',
                       midpoint=0.5, limits=c(0,1)) +
  scale_x_continuous(breaks=c(-73.15,-73.11),
                     labels=c(-73.15, -73.11)) +
  ylab("Longitude") + xlab("Latitude") +
  labs(fill=NULL) + ggtitle("2013")


postscript("Fig6.eps",height=8,width=11)
print(prob2013.map)
dev.off()


############################################################################################
## Figure 7: Map of selected census sector to display relative risk of dengue disease     ##
## for the period Jan 2009 - Dec 2015 (left panel), and specific temporal evolution of    ##
## the posterior mean estimates of relative risk and 95% credible intervals (right panel) ##
############################################################################################

## Right panel ##
mtl.IV <- data.frame(DengueData_INLA, Model.TLA.IV$summary.fitted.values)

selsec <- c("0527","1411","0102","1194","0320","0742","1238","0423")

mtl.IVshort <- mtl.IV[mtl.IV$sector %in% selsec,]
mtl.IVshort$sectorlabel <- mtl.IVshort$sector
mtl.IVshort$sectorlabel <- recode(mtl.IVshort$sectorlabel,
                                  "'0527'= c('Campo Hermoso');
                                  '1411'= c('Morrorrico');
                                  '0102'= c('Kennedy');
                                  '1194'= c('Provenza');
                                  '0320'= c('San Francisco');
                                  '0742'= c('Real de Minas');
                                  '1238'= c('Cabecera');
                                  '0423'= c('Girardot')")

theme_longitudinalRR <-	theme_bw() + theme(
  panel.border = element_rect(colour="black", size=0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=12),  
  axis.text.x  = element_text(size=10), 
  axis.title.y = element_text(size=12),
  axis.text.y  = element_text(size=12),
  plot.title   = element_text(size=36),
  legend.text  = element_text(size=12,lineheight=1.0),
  legend.title = element_text(size=12),
  legend.key.height=unit(1, "cm"),
  strip.text = element_text(face="bold", size=10),
  strip.background = element_rect(fill="white", colour="white")) 

RiskSelected <- ggplot(data=mtl.IVshort,aes(x=ID.year,y=mean))+
  geom_line() +
  geom_ribbon(aes(ymin=X0.025quant, ymax=X0.975quant), alpha=0.2)+
  scale_fill_gradient2(low='green4',mid='white',high='red4',
                       limits=c(0,5),midpoint=0) +
  facet_wrap(~sectorlabel,ncol=2)+
  ylab(expression(paste(r[it]))) + xlab("Epidemiological year")+
  geom_hline(yintercept = 1,linetype="dashed") + 
  scale_x_continuous(breaks=c(1,14,27,40,53,66,79),
                     labels=c("2009","2010","2011","2012","2013","2014","2015")) +
  theme_longitudinalRR


## Left panel ##
id.selected  <- unique(mtl.IVshort$sector)
id.selected.frame <- data.frame(id.selected = id.selected,mark=0)

Carto.sector.dataA <- merge(Carto.sector.np@data, id.selected.frame,
                            by.x="ID",by.y="id.selected",all.x=TRUE)
Carto.sector.dataA[is.na(Carto.sector.dataA)] <- 1
Carto.sector.dataA$id                         <- 1:94
Carto.sector.dataA$mark.factor                <- factor(Carto.sector.dataA$mark)

centroids.df       <- as.data.frame(coordinates(Carto.sector.np))
Carto.sector.dataB <- data.frame(Carto.sector.dataA,centroids.df)

id.selectedName <- data.frame(sector = id.selected)
id.selectedName$sectorName <- c('Morrorrico',
                                'Provenza',
                                'Real de Minas',
                                'San Francisco',
                                'Cabecera',
                                'Campo Hermoso',
                                'Girardot',
                                'Kennedy')

Carto.sector.dataC <- merge(Carto.sector.dataB,id.selectedName,
                            by.x="ID",by.y="sector",all.x=TRUE)

themeSelectedShort <-	theme_bw() + theme(
  panel.border = element_rect(colour="black", size=0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=14),  
  axis.text.x  = element_text( size=14.0), 
  axis.title.y = element_text(size=14),
  axis.ticks   = element_blank(),
  axis.text.y  = element_text( size=14),
  plot.title   = element_text(size=14),
  legend.text  = element_blank(),
  legend.title = element_blank(),
  strip.text   = element_text(face="bold", size=14),
  strip.background = element_rect(fill="white", colour="white")) 

base.line1 <- data.frame(x=c(-73.1390,-73.12) ,y=c(7.1529,7.169))
base.line2 <- data.frame(x=c(-73.1229,-73.16) ,y=c(7.1286,7.157))
base.line3 <- data.frame(x=c(-73.1347,-73.163),y=c(7.1217,7.145))
base.line4 <- data.frame(x=c(-73.1382,-73.158),y=c(7.1040,7.075))
base.line5 <- data.frame(x=c(-73.1246,-73.14),y=c(7.1022,7.07))
base.line6 <- data.frame(x=c(-73.1209,-73.10),y=c(7.0806,7.07))
base.line7 <- data.frame(x=c(-73.1084,-73.095),y=c(7.1170,7.084))
base.line8 <- data.frame(x=c(-73.0989,-73.105),y=c(7.1304,7.15))

base.text <- data.frame(
  x=c(-73.12,-73.161,-73.163,-73.158,-73.14,-73.10,-73.097,-73.105),
  y=c(7.172,7.16,7.148,7.074,7.068,7.068,7.082,7.155),
  sectorName=c("Kennedy","San Francisco","Girardot","Campo Hermoso",
               "Real de Minas","Provenza","Cabecera","Morrorrico"))

mapSelected <- ggplot() +
  geom_map(data=Carto.sector.dataB, map=map.sector.np, aes(fill=mark.factor, 
                                                           map_id=ID),color="black",size=0.08) +
  expand_limits(x=map.sector.np$long,y=map.sector.np$lat) +
  holes +
  themeSelectedShort + 
  geom_line(aes(x=x,y=y),data=base.line1,col="blue",size=1)+
  geom_line(aes(x=x,y=y),data=base.line2,col="blue",size=1)+
  geom_line(aes(x=x,y=y),data=base.line3,col="blue",size=1)+
  geom_line(aes(x=x,y=y),data=base.line4,col="blue",size=1)+
  geom_line(aes(x=x,y=y),data=base.line5,col="blue",size=1)+
  geom_line(aes(x=x,y=y),data=base.line6,col="blue",size=1)+
  geom_line(aes(x=x,y=y),data=base.line7,col="blue",size=1)+
  geom_line(aes(x=x,y=y),data=base.line8,col="blue",size=1)+
  geom_text(aes(x=x,y=y,label=sectorName), data = base.text, size=5,col="blue")+
  guides(fill=FALSE) +
  scale_fill_manual(values=c("pink","white")) +
  scale_x_continuous(breaks=c(-73.15,-73.11),
                     labels=c(-73.15, -73.11)) +
  ylab("Longitude") + xlab("Latitude") +
  labs(fill=NULL)  

cairo_ps("Fig7.eps",height=8,width=12)
grid.arrange(mapSelected,RiskSelected,nrow=1)
dev.off()