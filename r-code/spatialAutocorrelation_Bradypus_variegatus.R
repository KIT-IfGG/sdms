library(raster)
library(rgdal)
library(dismo)
library(rJava)
#library(maptools)
library(ncf) # correlogramm (there are alternatives!)
library(spdep) # neighbours
library(classInt)

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
backgr <- randomPoints(predictors, 1000)

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
plot(lat ~ lon, data=bradypus_pa, col=colorRampPalette(c("white", "red"))(nc)[ints], cex=pt.cex, pch=16)
points(lat ~ lon, data=bradypus_pa[bradypus_pa$pa == 1, ])

### For P data 
nc <- 20
brks <- classIntervals(bradypus_pa$residuals[bradypus_pa$pa == 1], nc, style = "equal")
ints <- findInterval(bradypus_pa$residuals[bradypus_pa$pa == 1], brks$brks)
plot(lat ~ lon, data=bradypus_pa[bradypus_pa$pa == 1,], col=colorRampPalette(c("white", "red"))(nc)[ints], cex=pt.cex, pch=16)

### Correlogram ####
cg <- correlog(bradypus_pa$lon, bradypus_pa$lat, bradypus_pa$residuals, increment=200, resamp=100, na.rm = TRUE, latlon = TRUE)
plot(cg$mean.of.class, cg$correlation, pch=16, col=ifelse(cg$p < 0.05, "red", "black"), xlab = "Distance", ylab="Spatial autocorrelation")
abline(h=0, col="blue")
legend("topright", legend=c("Sign.", "non sign."), col=c("red", "black"), pch=16)

