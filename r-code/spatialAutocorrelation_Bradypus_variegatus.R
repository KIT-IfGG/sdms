library(raster)
library(rgdal)
library(dismo)
library(rJava)
#library(maptools)
library(ncf) # correlogramm (there are alternatives!)
library(spdep) # neighbours
library(classInt)
library(corrplot)
library(nlme)

############################################### #
### Get data ####
############################################### #
### Path to the example species data.
file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
file
### Read the file
bradypus <- read.table(file, header=TRUE, sep=",")
head(bradypus)

### We only need columns 2 and 3:
bradypus <- bradypus[,2:3]
head(bradypus)

### Get environmental data  ####
data(wrld_simpl)   ### Data provided by "maptools". Same as by getData, GADM.

### Environmental data, predictirs ####
### BioClim data provided with the package dismo.
files <- list.files(path=paste(system.file(package="dismo"),"/ex", sep=""), pattern="grd", full.names=TRUE )
predictors <- stack(files)

############################################### #
### Model ####
############################################### #

### Create background points
backgr <- randomPoints(predictors, 500)

### MaxEnt with PB-data (often also refered to as presence-only data)
me <- maxent(predictors[[c("bio1", "bio12")]], bradypus, backgr)

### Create dataset for residual analysis ####
colnames(backgr) <- colnames(bradypus)
bradypus_pa <- rbind(bradypus, backgr)
bradypus_pa <- data.frame(bradypus_pa, pa = c(rep(1, nrow(bradypus)), rep(0, nrow(backgr))))
bios_pa <- extract(predictors, bradypus_pa[,1:2])
bradypus_pa$fitted <- predict(me, bios_pa)

bradypus_pa$residuals <- bradypus_pa$pa  - bradypus_pa$fitted 

boxplot(fitted ~ pa, bradypus_pa)

### Plot residuals in space ####
### For all data
nc <- 10
pt.cex <- 0.5
brks <- classIntervals(bradypus_pa$residuals, nc, style = "quantile")
ints <- findInterval(bradypus_pa$residuals, brks$brks)
plot(lat ~ lon, data=bradypus_pa, col=colorRampPalette(c("blue", "red"))(nc)[ints], cex=pt.cex, pch=16)
points(lat ~ lon, data=bradypus_pa[bradypus_pa$pa == 1, ])

### Correlogram ####
cg <- correlog(bradypus_pa$lon, bradypus_pa$lat, bradypus_pa$residuals, increment = 200, resamp = 100, na.rm = TRUE, latlon = TRUE)
plot(cg$mean.of.class, cg$correlation, pch=16, col=ifelse(cg$p < 0.01, "red", "black"), xlab = "Distance", ylab="Spatial autocorrelation")
abline(h=0, col="blue")
legend("topright", legend=c("Sign.", "non sign."), col=c("red", "black"), pch=16)
cg$n

### Account for spatial autocorrelation ####
### Use GLS  and correlation structure ####
dat <- rbind(bradypus, backgr)
dat$pa <- c(rep(1, nrow(bradypus)), rep(0, nrow(backgr)))
biodat <- extract(predictors[[c("bio1", "bio12")]], dat[,1:2])
dat <- cbind.data.frame(dat, biodat)

### Spatial project ####
xy <- SpatialPoints(dat[,c("lon", "lat")], CRS(projection(predictors)))
xy <- spTransform(xy, CRS("+init=epsg:2317"))
xy <- spTransform(xy, CRS("+init=epsg:6610"))

plot(xy)

predictors_m <- projectRaster(predictors, crs=CRS("+init=epsg:6610"))

dat$x <- xy@coords[,1]
dat$y <- xy@coords[,2]

m1 <- gls(pa ~ poly(bio1, 2) * poly(bio12, 2), data = dat)
m2 <- gls(pa ~ poly(bio1, 2) * poly(bio12, 2), data = dat, correlation = corGaus(form = ~ x + y))
m3 <- gls(pa ~ poly(bio1, 2) * poly(bio12, 2), data = dat, correlation = corLin(form = ~ x + y))

### The model with the lowest AIC or BIC is preferred ###
AIC(m1)
AIC(m2)
AIC(m3)

BIC(m1)
BIC(m2)
BIC(m3)   

rbind(coef(m1), coef(m2), coef(m3))

round(cor(cbind(pa=dat$pa, m1=fitted(m1), m2=fitted(m2), m3=fitted(m3)), method = "spearman"), 2)
corrplot(cor(cbind(pa=dat$pa, m1=fitted(m1), m2=fitted(m2), m3=fitted(m3)), method = "spearman"), type = "lower")
### -> Best BIC does not have best correlation!

### m2 ###
summary(m2)  ### "range" - in m
corr_range <- 16725.6
dat$residuals <- resid(m2)

ncf.cor <-  correlog(dat$x, dat$y, dat$residuals, increment = 2000, resamp = 100, na.rm = TRUE, latlon = FALSE)
#ncf.cor$n
plot(ncf.cor$mean.of.class, ncf.cor$correlation, type = "b", pch = 16, col=ifelse(ncf.cor$p < 0.01, "red", "black"), xlab = "Distance", ylab = "Spatial autocorrelation", xlim=c(0, 400000))
abline(h = 0, col = "blue")
abline(v = corr_range , lty = 2)
legend("bottomright", legend=c("Sign.", "non sign."), col=c("red", "black"), pch=16)

### m3 ####
summary(m3)
corr_range_lin <- 2856692

dat$residuals <- resid(m3)
ncf.cor <- correlog(dat$x, dat$y, dat$residuals, increment = 2000, resamp = 100, na.rm = TRUE, latlon = FALSE)
#ncf.cor$n
plot(ncf.cor$mean.of.class, ncf.cor$correlation, type = "b", pch = 16, col=ifelse(ncf.cor$p < 0.01, "red", "black"), xlab = "Distance", ylab = "Spatial autocorrelation", xlim=c(0, 3000000))
abline(h = 0, col = "blue")
abline(v = corr_range_lin , lty = 2)
legend("bottomright", legend=c("Sign.", "non sign."), col=c("red", "black"), pch=16)

### Maps ####
e1 <- evaluate(dat[dat$pa == 1,c("lon", "lat")], dat[dat$pa == 0,c(c("lon", "lat"))], m1, predictors)
e3 <- evaluate(dat[dat$pa == 1,c("lon", "lat")], dat[dat$pa == 0,c(c("lon", "lat"))], m3, predictors)

thr1 <- threshold(e1)$kappa
thr3 <- threshold(e3)$kappa

pred1 <- predict(predictors_m[[c("bio1", "bio12")]], m1)
pred3 <- predict(predictors_m[[c("bio1", "bio12")]], m3)
#pred_me <- predict(predictors_m[[c("bio1", "bio12")]], me)


pred1[pred1[] < thr1] <- NA
pred3[pred3[] < thr3] <- NA

pred1[] <- (pred1[] - min(pred1[], na.rm = TRUE)) / (max(pred1[], na.rm = TRUE) - min(pred1[], na.rm = TRUE))
pred3[] <- (pred3[] - min(pred3[], na.rm = TRUE)) / (max(pred3[], na.rm = TRUE) - min(pred3[], na.rm = TRUE))

x11(width =9, height = 6)
par(mfrow=c(1,2))
plot(pred1, col = colorRampPalette(c("green", "darkgreen"))(10))
plot(pred3, col = colorRampPalette(c("green", "darkgreen"))(10))
