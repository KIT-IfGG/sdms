### Author: Klara Dolos
### Calculate species distribution models using R package dismo.
#install.packages(c("raster", "rgdal", "dismo", "rJava", "gbm", "maptools"))
### To run MaxEnt you need to install the packages below. You might need to
### Install or update also Java. You need to READ THE ERROR MESSAGES R provides you
### and follow the suggestions made. You will probably need some time to get
### through all of this, so do this at home with a good internet connection!
### When these packages can be loaded, execute the code untill the 
### flag "TESTING DONE" without errors.
library(raster)
library(rgdal)
library(dismo)
library(rJava)
library(gbm)
library(maptools)

library(corrplot)
library(hier.part)


### Part I ####
### Get species data ####
### Path to the example data.
file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
file
# read it
bradypus <- read.table(file, header=TRUE, sep=",")
head(bradypus)

# we only need columns 2 and 3:
bradypus <- bradypus[,2:3]
head(bradypus)


### Get environmental data  ####
data(wrld_simpl)   ### Data provided by "maptools". Same as by getData, GADM.
x11()
plot(wrld_simpl)

### Environmental data, predictors ####
### BioClim data provided with the package dismo.
files <- list.files(path=paste(system.file(package="dismo"),"/ex", sep=""), pattern="grd", full.names=TRUE )
predictors <- stack(files)
plot(predictors)

plot(predictors, 1)
plot(wrld_simpl, add=TRUE)
points(bradypus, col="blue")

### TESTING DONE ####

### Create background points
set.seed(8)
backgr <- randomPoints(predictors, 500)

x11()
plot(wrld_simpl)
points(bradypus, col="blue", cex=0.3, pch=16)
points(backgr, col="red", cex=0.3, pch=16)

### Extracting values from rasters
presvals <- extract(predictors, bradypus)
absvals <- extract(predictors, backgr)
head(presvals)

### Create predictor dataset
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
sdmdata[,"biome"] <- as.factor(sdmdata[,"biome"])
head(sdmdata)
tail(sdmdata)

### Correlation among predictors
x11()
pairs(sdmdata[,2:8], cex=0.1, fig=TRUE)
corrplot(cor(sdmdata[,1:8]), type= "lower", diag=FALSE)

### Model fitting using different model types ####
### Algorithms require different data structure, i.e. sdmdata, predictors, bradypus.

### Logistic regression (should be used with PA-data!)
m1 <- glm(pb ~ bio1 + bio5 + bio12, data=sdmdata, family=binomial)
summary(m1)

m2 <- step(glm(pb ~ ., data=sdmdata))
summary(m2)

### Climatic envelope model with presence-only data (P data)
bc <- bioclim(presvals[,c("bio1", "bio5", "bio12")])
bc

### MaxEnt with PB-data (often also refered to as presence-only data)
me <- maxent(predictors, bradypus)
me

### Boosted regression tree with pa-data
brt <- gbm(pb ~ bio1 + bio5 + bio12, data=sdmdata, distribution="bernoulli")
brt

### Prediction ####
p.m1 <- predict(predictors, m1, type="response")
p.bc <- predict(predictors, bc)
p.me <- predict(predictors, me)
p.brt <- predict(predictors, brt, n.trees=10, type="response")

x11()
par(mfrow=c(2,2))
plot(p.m1, main=class(m1)[1])
plot(p.bc, main=class(bc)[1])
plot(p.me, main=class(me)[1])
plot(p.brt, main=class(brt)[1])
dev.off()

graphics.off()  ### Close all graphic devices

### Part II ####

### Model evaluation and quality ####
ths <- seq(0,1, len=50)
e <- evaluate(p=presvals, a=absvals, bc, tr=ths)

x11()
par(mfrow=c(1,3))
plot(e,'ROC')
boxplot(e, col=c("blue", "red"), main="Probs")
density(e)
dev.off()

### Exkurs
dat <- sdmdata
dat$pred <- predict(bc, dat)   ### Calculate modelled "probability" of occurrance. "Prediction".
names(dat)
dat[1:10,c("pb", "pred")]

x11(width=8, height=6)
par(mfrow=c(1,2))
boxplot(pred ~ pb, dat, col=c("red", "blue"))
grid()
#stripchart(pred ~ pb, dat)
pcol <- rgb(0,0,1,0.4)
acol <- rgb(1,0,0,0.6)
hist(dat$pred[dat$pb==0], col=acol, border="red", freq=F, main="")
hist(dat$pred[dat$pb==1], add=TRUE, , col=pcol, border="white" freq=F)
legend("topright", legend=c("Presences", "Absences"), col=c(pcol, acol), pch=15, pt.cex=2)
grid()
box()
dev.off()

### Similar to the density plot. Distribution of predicted occurrence probabilities.
x11()
par(mfrow=c(1,2))
hist(dat$pred[dat$pb ==1])
hist(dat$pred[dat$pb ==0])

### End Exkurs

### Add labels to find threshold
str(e)
hlpr <- cbind.data.frame(FPR=e@FPR, TPR=e@TPR,  labels=ths)
hlpr <- round(hlpr, 2)
d <- duplicated(hlpr[,c("FPR", "TPR")])
table(d)
hlpr <- hlpr[!d,]

x11()
par(mfrow=c(1,1))
plot(e,'ROC', type="l")
text(hlpr$FPR, hlpr$TPR,  labels=hlpr$labels)

th <- 0.3   ### Not a nice example...

### Automatic way to find the best threshold
ths <- threshold(e); ths
th <- ths$kappa   ### e.g. use kappa statistics

### Plot species distribution using a threshold ####
mycolors <- c("lightgrey", "blue")
bradypus_distribution <- p.bc 
bradypus_distribution[bradypus_distribution[]<th] <- 0
bradypus_distribution[bradypus_distribution[]>=th] <- 1
x11()
plot(bradypus_distribution, main=paste0("Threshold = ", round(th,2)), col=mycolors, legend=FALSE)
legend("topright", legend=c("Bradypus habitat", "No habitat"), col=mycolors, pch=15)

### Variable importance ####

### GLM: Numbers given in the table and figure.
hp <- hier.part(sdmdata$pb, sdmdata[,c("bio1", "bio5", "bio12")], family=binomial)
box()
hp$I.perc

### MaxEnt: "Analysis of variable contributions"; in your browser.
me

### BRT: Numbers given in the table and figure.
summary(brt) 
str(summary(brt))

### BioClim: ...?
summary(bc) 

### Part III ####

### Validation: k-fold data partitioning for BC ####
pres <- sdmdata[sdmdata[,1] == 1, 2:9]
back <- sdmdata[sdmdata[,1] == 0, 2:9]

k <- 5
group <- kfold(pres, k)

e <- list()
for (i in 1:k) {
  train <- pres[group != i,]    # calibration dataset, training data
  test <- pres[group == i,]    # test dataset, validation dataset
  bc <- bioclim(train)
  e[[i]] <- evaluate(p=test, a=back, bc)
}

auc <- sapply(e, function(x) x@auc)
round(mean(auc),2)
round(sd(auc),2)

round(median(auc),2)
round(cv(auc),2)

### Validation: k-fold data partitioning for ME incl. preditions ####
bradypus   ### presences
background <- randomPoints(predictors[[1]], 1000) ### background ("absences")
predictors ### "Bioclim dataset"

k <- 5
group <- kfold(bradypus, k)

e <- list()
for (i in 1:k) {
  train <- bradypus[group != i,]    # calibration dataset, training data
  test <- bradypus[group == i,]    # test dataset, validation dataset
  me <- maxent(predictors, train, background)
  e[[i]] <- evaluate(p = test, a = background, model = me, x = predictors)
  
  p.me <- predict(predictors, me)
  writeRaster(p.me, file = paste0("figures/validation_me/pred_", i, ".tif"), overwrite=TRUE)
}

auc <- sapply(e, function(x) x@auc)
round(mean(auc),2)
round(sd(auc),2)

round(median(auc),2)
round(cv(auc),2)

### Prediction uncertainty based on k-fold cross-validation ####
rfis <- list.files("figures/validation_me", full = T)
p.mes <- stack(rfis)
cv.map <- calc(p.mes, cv)
plot(cv.map, col=heat.colors(5), main="CV for predictions")
