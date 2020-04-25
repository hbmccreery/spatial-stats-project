library(fields)
library(geoR)

## MERCURY DATA
# read data
AMNet.ALL <- read.csv("AMNet-ALL.csv", stringsAsFactors=FALSE)
SiteID_info <- read.csv("SiteID_info.csv", stringsAsFactors=FALSE)

# select data within  year (2016)
AMNet.ALL$CollEnd <- as.Date(AMNet.ALL$CollEnd)
AMNet.ALL$CollStart <- as.Date(AMNet.ALL$CollStart)
oneYear <- AMNet.ALL[AMNet.ALL$CollStart >= as.Date('2015-01-01') & AMNet.ALL$CollEnd < as.Date('2017-01-01'),]

# select observations with "valid" data as defined in metadata
validData <- oneYear[(oneYear$PBMVal == 'A' | oneYear$PBMVal == 'B') &
                       (oneYear$GOMVal == 'A' | oneYear$GOMVal == 'B') &
                       (oneYear$GEMVal == 'A' | oneYear$GEMVal == 'B'),]

# drop missing/invalid data (represented by -9)
validData <- validData[validData$PBM != -9 & validData$GOM != -9 & validData$GEM != -9,]

# find means over this data
# tried to use cast/reshape but threw weird errors I couldn't debug
stations <- unique(validData$SiteID)
numStations = length(stations)
means <- data.frame(Station = character(numStations),
                    PBM = double(numStations),
                    GOM = double(numStations),
                    GEM = double(numStations),
                    stringsAsFactors=FALSE)

i <- 1
for (station in stations){
  stationData <- validData[validData$SiteID == station,]
  means[i,] <- c(station, mean(stationData$PBM), mean(stationData$GOM), mean(stationData$GEM))
  i <- i + 1
}

# clean names of things
colnames(means) <- c('Site.ID', 'PBM', 'GOM', 'GEM')
means$PBM <- as.double(means$PBM)
means$GOM <- as.double(means$GOM)
means$GEM <- as.double(means$GEM)

# create df with site location info, drop hawaii b/c too far away
dev.off()
withLoc <- merge(means, SiteID_info)
withLoc <- withLoc[withLoc$Site.ID != 'HI00',]

# PLOTTING
# pbm
quilt.plot(withLoc[,c(9,8)], withLoc$PBM, nx=40, ny=40)
US(add=TRUE)
mtext('Average Particulate Bound Mercury (' ~ pg/m^3 ~ ') in 2016')
dev.print(png, '~/R/PBM.png', width=800)

# gom
quilt.plot(withLoc[,c(9,8)], log(withLoc$GOM), nx=40, ny=40)
US(add=TRUE)
mtext('Log scaled Average Gaseous Oxidized Mercury (' ~ pg/m^3 ~ ') in 2016')
dev.print(png, '~/R/GOM.png', width=800)

# gem
quilt.plot(withLoc[,c(9,8)], withLoc$GEM, nx=40, ny=40)
US(add=TRUE)
mtext('Average Gaseous Elemental Mercury (' ~ ng/m^3 ~ ') in 2016')
dev.print(png, '~/R/GEM.png', width=800)

## PUBLIC HEALTH DATA
# for plotting later
library(usmap)
library(ggplot2)

# read in pre-processed data
lungCounty = read.csv('~/R/lungCounty.csv')

plotting <- map_with_data(lungCounty, values='AGE_ADJUSTED_RATE')
plotting$AGE_ADJUSTED_RATE <- as.double(plotting$AGE_ADJUSTED_RATE)

plot_usmap(data=plotting, values='AGE_ADJUSTED_RATE') + scale_fill_continuous(name = "Lung Cancer Incidence Rate", label = scales::comma) + theme(legend.position = "right")  + labs(title = "Lung Cancer Incidence Rate in the US")
dev.print(png, '~/R/GEM.png', width=800)

#VARIOGRAMS
GEM.variog <- variog(coords = withLoc[,c(9,8)], data=withLoc$GEM, estimator.type = 'classical')
GOM.variog <- variog(coords = withLoc[,c(9,8)], data=withLoc$GOM, estimator.type = 'classical')
PBM.variog <- variog(coords = withLoc[,c(9,8)], data=withLoc$PBM, estimator.type = 'classical')

plot(GEM.variog, main='GEM Variogram')
lines(GEM.variog)
plot(GOM.variog, main='GOM Variogram')
lines(GOM.variog)
plot(PBM.variog, main='PBM Variogram')
lines(PBM.variog)

## MLE ESTIMATOR
filtered <- withLoc[withLoc$Longitude > -100,]
dist.mat <- rdist.earth(filtered[,c(9,8)])
z <- log(filtered$GOM) - mean(log(filtered$GOM))

# grd.size <- 30
# rand.grd.long <- seq(min(filtered$Longitude), max(filtered$Longitude), length.out=grd.size)
# rand.grd.lat <- seq(min(filtered$Latitude), max(filtered$Latitude), length.out=grd.size)
# rand.grd <- as.matrix(expand.grid(rand.grd.long,rand.grd.lat))
# dist.mat <- rdist.earth(rand.grd)
# 
# Sigma <- Matern(dist.mat,range=0.2,nu=1.5)
# set.seed(10)
# z <- t(chol(Sigma)) %*% rnorm(grd.size^2)
# 
# quilt.plot(rand.grd, z, nx=grd.size, ny=grd.size)
# US(add=TRUE)

ell1 <- function(p){
  # p = (tau, sigma^2, a)
  Sigma.c <- t(chol(diag(exp(p[1])^2,length(z)) + Matern(dist.mat,range=p[3],nu=0.5)*p[2])) # lower triangular
  out <- forwardsolve(Sigma.c,z)
  quad.form <- sum(out^2)
  det.part <- 2*sum(log(diag(Sigma.c)))
  
  # up to the normalizing constants
  0.5*det.part + 0.5*quad.form
}

# ell2 <- function(p){
#   # as a function of a,nu
#   Sigma.c <- t(chol(Matern(dist.mat,range=p[1],nu=0.5)*p[2])) # lower triangular
#   out <- forwardsolve(Sigma.c,z)
#   quad.form <- sum(out^2)
#   det.part <- 2*sum(log(diag(Sigma.c)))
#   
#   # up to the normalizing constants
#   0.5*det.part + 0.5*quad.form
# }

param.grd.size <- 50
a_min <- 0.01
a_max <- 1000
var_min <- 1
var_max <- 5
tau_min <- 0
tau_max <- 1

a <- seq(a_min,a_max,length.out=param.grd.size)
var_seq <- seq(var_min,var_max,length.out=param.grd.size)
tau_seq <- seq(tau_min,tau_max,length.out=param.grd.size)

par.grd <- as.matrix(expand.grid(tau_seq,var_seq,a))

ells <- NULL
for(i in 1:dim(par.grd)[1]){
  ells[i] <- ell1(c(par.grd[i,1],par.grd[i,2],par.grd[i,3]))
  print(i/dim(par.grd)[1])
}

## 2d grid search plot
par(mfrow=c(1,1))
image.plot(a,var_seq,matrix(ells,param.grd.size,param.grd.size),xlab="range",ylab="var")
this <- order(ells)[1]
points(rbind(par.grd[this,]),pch=19,col="white")

## 3d grid search plot
max_val <- which.min(ells)
par.grd[max_val,]
these <- par.grd[,1] == par.grd[max_val,1]
subset_grd <- par.grd[these,c(2,3)]
subset_ells <- ells[these]
quilt.plot(subset_grd, subset_ells, nx=param.grd.size, ny=param.grd.size, xlab='var', ylab='lambda')

this <- order(subset_ells)[1]
points(rbind(subset_grd[this,]),pch=19,col="white")

out <- optim(par=c(0,1,500),fn=ell1,hessian=TRUE,method="L-BFGS-B",lower=c(tau_min, var_min, a_min),
             upper=c(tau_max, var_max, a_max))

## MLEs
MLEs <- out$par
tau <- MLEs[1]
sigma2 <- MLEs[2]
a <- MLEs[3]

MLEs

## MLEs with standard errors
I.inv <- solve(out$hessian)
SEs <- sqrt(diag(I.inv))

CIs <- cbind(out$par - 1.96*SEs,out$par + 1.96*SEs)
rownames(CIs) <- c("tau","sigma2","a")
colnames(CIs) <- c("2.5 %","97.5 %")
CIs # truth is (0.2, 1.5)

##### KRIGING TO A GRID
## Create grid to predict on
lonvals <- seq(min(filtered$Longitude), max(filtered$Longitude), length.out = 100)
latvals <- seq(min(filtered$Latitude), max(filtered$Latitude), length.out = 250)
pred.grd <- as.matrix(expand.grid(lonvals, latvals))

## Matrices for kriging
dist0.mat <- rdist.earth(filtered[,c(9,8)], pred.grd)
Sigma0 <- (sigma2 * Matern(dist0.mat, range=a, nu=0.5))
#diag(exp(tau)^2,length(pred.grd)) + 
dist.mat <- rdist.earth(filtered[,c(9,8)])
Sigma <- sigma2 * Matern(dist.mat, range=a, nu=0.5)

## Simple kriging predictor
sk <- t(Sigma0) %*% solve(Sigma) %*% z

par(mfrow=c(1,2))
quilt.plot(filtered[,c(9,8)],z + mean(log(filtered$GOM)))
US(add=TRUE)
quilt.plot(pred.grd,sk + mean(log(filtered$GOM)))
US(add=TRUE)

## Uncertainties
sk.var <- sigma2 - diag(t(Sigma0) %*% solve(Sigma) %*% Sigma0)
sk.se <- sqrt(sk.var)
par(mfrow=c(1,1))
quilt.plot(pred.grd, sk.se)
US(add=TRUE)
