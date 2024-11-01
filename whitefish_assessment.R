## Development branch of SPiCT needed (spinup option)
##remotes::install_github("DTUAqua/spict/spict", ref = "dev")
library(spict)

##surv = read.csv("Survey_mesh40_44.csv", sep=";", dec=",")
surv = read.csv("Survey.csv", sep=";", dec=",")


catch = read.csv("catches.csv",sep=";",dec=",")

##effort = read.csv("Effort.csv", sep=";", dec=".")
## This yearly effort is not accurate - seems like different catchability between quarters
## and in Q4 effort is only partially recorded so should not be used!


Catch_quarter = read.csv2("Catches_quarterly.csv")

Catch_quarter$Value <- as.numeric(Catch_quarter$Value)
Catch_quarter <- Catch_quarter[Catch_quarter$Quarter!= "total",]

Effort_quarter <- read.csv("Effort_quarterly.csv", sep=";")

Catch_quarter$time <- as.numeric(Catch_quarter$Year)+as.numeric(Catch_quarter$Quarter)/4-0.25
Catch_quarter_u <- Catch_quarter[order(Catch_quarter$time),] 

Effort_quarter$time <- as.numeric(Effort_quarter$Year)+as.numeric(Effort_quarter$Quarter)/4-0.25
Effort_quarter_u <- Effort_quarter[order(Effort_quarter$time),] 

sel = which(Catch_quarter_u$JYear >= min(Effort_quarter$Year) )
selC = which(catch[,1]>=1997 & catch[,1]<=2022)  

dataq <- data.frame(Year = Catch_quarter_u$Year[sel],
                    cpue = Catch_quarter_u$Value[sel] / Effort_quarter_u$Value*1e3,
                    time = Catch_quarter_u$time[sel],
                    catch =Catch_quarter_u$Value[sel],
                    effort = Effort_quarter_u$Value,
                    Quarter = Catch_quarter_u$Quarter[sel]
                    )

dataq$effort[ dataq$Quarter==4 ] = NA
dataq$cpue[ dataq$Quarter==4 ] = NA

datQ <- split(dataq,dataq$Quarter)

inp = list( obsI = list(surv[,5],datQ[[2]]$cpue,datQ[[3]]$cpue), timeI = list(surv[,1]+0.75,datQ[[2]]$time+1/8,datQ[[3]]$time+1/8), obsC = catch[selC,2], timeC=catch[selC,1])

## Use survey + Q2 and Q3 commercial CPUE.
## Q1 cpue has some outliers and most catches are taken in Q2 and Q3,
## which is is why only these are included.

inps <- list()
models <- list()

inps[[1]] <- inp

## Prior means and sds
## r is often estimated unrealisticly high for this dataset
## Fishbase says "medium" resilience, i.e. r between 0.2 and 1.0
## https://www.fishbase.se/summary/Coregonus-wartmanni.html
logrmean <- log(0.4)
logrsd <- 0.4
logsdbmean <- log(0.15)
logsdbsd <- 0.4
logsdcmean <- log(0.1)
logsdcsd <- 0.2


rrange <- exp(c( logrmean-2*logrsd, logrmean+2*logrsd))


modelNr <- 1
## fix to schaefer n=2
inps[[modelNr]]$phases$logn=-1
inps[[modelNr]]$ini$logn=log(2)

inps[[modelNr]]$priors$logr <- c(logrmean,logrsd,1)  

## disable default priors
inps[[modelNr]]$priors$logalpha=c(0,0,0)
inps[[modelNr]]$priors$logbeta=c(0,0,0)

models[[modelNr]] = fit.spict(inps[[modelNr]])



modelNr <- 2
## As 1, but include spin-up
## Spin-up means starting the model earlier (here 10 years)
## than the data - this ensures that the first biomass state (after removing spin-up)
## cannot be way above K


inps[[modelNr]] <- inps[[1]]
inps[[modelNr]]$nspinup = 160
models[[modelNr]] = fit.spict(inps[[modelNr]])

#### Add sdb and sdc priors (drop this as it is harder to argue for these priors in the paper).
##modelNr <- 3
##inps[[modelNr]] = inps[[2]]
##inps[[modelNr]]$priors$logsdb <- c(logsdbmean,logsdbsd,1)
##inps[[modelNr]]$priors$logsdc=c(logsdcmean,logsdcsd,1)

##models[[modelNr]] = fit.spict(inps[[modelNr]])


## Fix to (almost) fox
modelNr <- 3
inps[[modelNr]] = inps[[2]]
inps[[modelNr]]$ini$logn = log(1.1) ## almost fox
inps[[modelNr]]$priors$logr <- c(log(exp(logrmean)*1.1/2),logrsd,1) ## adjust r prior (reference: Kokkalis et. al 2024: Good practices for surplus production models)

models[[modelNr]] = fit.spict(inps[[modelNr]])

## Residuals seem to increase over time, try increase uncertainty of low observations
modelNr <- 4
inps[[modelNr]] = inps[[2]]
pp=0.1
for(ii in 1:length(inps[[modelNr]]$obsI)){
    inps[[modelNr]]$stdevfacI[[ii]] = 1/inps[[modelNr]]$obsI[[ii]]^pp
    inps[[modelNr]]$stdevfacI[[ii]] = inps[[modelNr]]$stdevfacI[[ii]] / mean(inps[[modelNr]]$stdevfacI[[ii]])
}
models[[modelNr]] <- fit.spict(inps[[modelNr]])

## Try with shorter time-series
modelNr <- 5
inps[[modelNr]] = inps[[2]]
inps[[modelNr]]$nspinup=0
inps[[modelNr]] <- shorten.inp(inps[[modelNr]],mintime=2001,maxtime=Inf)
models[[modelNr]] = fit.spict(inps[[modelNr]])



## Add diagnostics (takes some time)
for(ii in 1:length(models)){
    models[[ii]] <- calc.osa.resid(models[[ii]])
    models[[ii]] <- retro(models[[ii]])
    models[[ii]] <- hindcast(models[[ii]])
}

mohnsrho <- list()
masevals <- list()
for(ii in 1:length(models)){
    mohnsrho[[ii]] <-  mohns_rho(models[[ii]])
    hcInfo <- spict:::extract.hindcast.info(models[[ii]])
    masevals[[ii]] <- hcInfo$mase[,2]
}

masevals
## all below 1, which is good
mohnsrho
plotspict.compare(models,plot.unc=FALSE)
## model 2 and 4 have least retro,
## and is somewhere "in the middle",
## so if we we want to present just one
## I would choose model 2.


plot(models[[2]])
plotspict.retro(models[[2]])
plotspict.diagnostic(models[[2]])

