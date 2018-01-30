### Install the R-package ###
#install.packages("rgbif")
### This might result in an error message. Install the suggested packages AND update other missing programms on you computer.
library(rgbif)
library(raster)
library(dismo)
library(maptools)
library(colorRamps)
library(classInt)

### Get species data ####
key <- name_suggest(q='Cicendia filiformis', rank='species')$key[1]
n <- occ_count(taxonKey=key, georeferenced=TRUE); n
presences <- occ_search(taxonKey=key, limit=n)
head(presences$data)

presences <- as.data.frame(presences$data[, c("decimalLatitude", "decimalLongitude")])
presences <- na.omit(presences)

### Save the data you will work with! Load from file from than on ####
saveRDS(presences, "data/cicendia_filiformis.rds")

### You better store this data part in a seperate R-Skript file named e.g. "get_my_species_data.R".

### Read in the data ###
presences <- readRDS("data/cicendia_filiformis.rds")

### Bioclim data at coarse resolution can be dowloaded within R
bioclim <- getData('worldclim', download=FALSE, path="data", var='bio', res=5) ### Use 'download = TRUE' at the first time! 
projection(bioclim)

data(wrld_simpl)  ### ?wrld_simpl
projection(wrld_simpl)

### Same projection but different projection string, fix
projection(wrld_simpl) <- projection(bioclim)

### ALWAYS plot the data to see if everything is ok!
plot(bioclim[[3]])

### Create spatial points dataframe
presences <- SpatialPoints(presences[,c("decimalLongitude", "decimalLatitude")], proj4string=CRS(projection(bioclim)))

x11()
plot(wrld_simpl)
points(decimalLatitude ~ decimalLongitude, data=presences, pch=16, cex=0.5, col="red")

### Use only europe; decide for your project, what to to!
europe <- wrld_simpl[wrld_simpl$REGION==150,]
europe@data 
europe <- europe[europe$NAME!="Russia",]

x11()
plot(europe)
points(decimalLatitude ~ decimalLongitude, data=presences, pch=16, cex=0.5, col="red")

saveRDS(europe, "data/euro_shp.rds")

### Exclude presence data outside the region ####
presences_europe <- presences[is.na(over(europe, presences)),]  

x11()
plot(europe)
points(presences_europe, pch=16, cex=0.5, col="red")

### Crop predictor raster dataset to study region ####
bioclim_europe <- crop(bioclim, europe)
bioclim_europe <- mask(bioclim_europe, europe)

saveRDS(bioclim_europe, "data/bioclim_euro.rds")

### Create background points for the model ####
background <- randomPoints(bioclim_europe, 20000)
#me <- maxent(bioclim_europe, presences_europe@coords, background)

### VARIABLE SELECTION ####
### Backward  selection +  ecological knowledge from literature. It is often resonable to use a temperature related and a precipitation related variable and not only t or p variables.
### Use variable importance and AUC for variable selection. Account for multicollinearity!

### You might like to try this:
#install.packages("MaxentVariableSelection")
#library(MaxentVariableSelection)
#vignette("MaxentVariableSelection")

me1 <- maxent(bioclim_europe, presences_europe@coords, background)
me <- maxent(bioclim_europe[[c("bio10", "bio18")]], p=presences_europe@coords, a=background) 

### Arguments of maxent: https://groups.google.com/forum/#!topic/maxent/yRBlvZ1_9rQ

### Model evaluation ####

par(mfrow=c(1,2))
plot(me1)
plot(me)

e1 <- evaluate(presences_europe, background, me1, bioclim_europe)
e <- evaluate(presences_europe, background, me, bioclim_europe)
e1
e

### Let's assume we have done a valid variable selection and evaluate the final model only
thr <- threshold(e)

### ROC-curve and density plot
par(mfrow=c(1,2)) 
plot(e, "ROC")
density(e)

### Response curves ####
response(me)

### 2D
np <- 30
newdata <- expand.grid(bio10=seq(145, 200, len=np), bio18=seq(0, 240, len=np))
newdata$pred <- predict(me, newdata)

### Use threshold to show distribution
newdata$pred[newdata$pred<thr$sensitivity] <- NA

### Create classes of site suitability
cInt <- classIntervals((newdata$pred))

xdiff <-diff(unique(newdata$bio10))[1]
ydiff <-diff(unique(newdata$bio18))[1]

mypalette <- colorRampPalette(c("lightgreen", "darkgreen"))
newdata$colors <- findColours(cInt, mypalette(length(cInt$brks)))

par(mfrow=c(1,1), mar=c(5,5,1,1))
symbols(x=newdata$bio10, y=newdata$bio18, rectangles=matrix(rep(c(xdiff, ydiff), nrow(newdata)), ncol=2, byrow=T), bg=newdata$colors, fg="white", inches=F, xlab="Temperature of warmest quarter (°dC)", ylab="Precipitation of warmest quarter (mm)")
contour(x=unique(newdata$bio10), y=unique(newdata$bio18), z=matrix(newdata$pred, nrow=np), add=T, levels=unique(round(cInt$brks,1)), labcex = 1.3)
mtext("Cicendia filiformis", side=3, line=-1.3, font=3)
mtext(paste0("AUC = " , round(e@auc, 2), " "), side=1, line=-2.3, adj=1)
mtext(paste0("Pearson r² = " , round(e@cor, 2), " "), side=1, line=-1.3, adj=1)

### Q: Is this response curve reasonable? Why?

### Plot distribution map ####
pred <- predict(me, bioclim_europe)
plot(pred)
distr <- pred
distr[distr < thr$sensitivity] <- NA
# cInt <- classIntervals((newdata$pred))

plot(distr, col=mypalette(10), breaks=cInt$brks, legend=F)
points(presences_europe, pch=16, cex=0.1, col="black")
plot(europe, add=T)
mtext("Cicendia filiformis", side=3, line=-1.3, font=3)
mtext(paste0("AUC = " , round(e@auc, 2), " "), side=1, line=-2.3, adj=1)
mtext(paste0("Pearson r² = " , round(e@cor, 2), " "), side=1, line=-1.3, adj=1)
### Add legend!

### Climate change projection ####
cc <- getData('CMIP5', var="bio", res=5, rcp=85, model='HD', year=70, download=TRUE, path="data")
cc <- crop(cc, bioclim_europe)
cc <- mask(cc, bioclim_europe)
names(cc) <- names(bioclim_europe)

pred_cc <- predict(me, cc)

distr_cc <- pred_cc
distr_cc[distr_cc[] < thr$sensitivity] <- NA

### Create figure for distribution for current and climate change projection ####
#x11(width=12)
pdf("figures/maps.pdf", width=12, pointsize = 16)
par(mfrow=c(1,2))
plot(distr, col=mypalette(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(europe, add=T)
mtext("Current climate ", 3,-2.2)
mtext("Cicendia filiformis ", 3,-1.2, font=3)
mtext("WGS84 ", 1,-1.2, adj=1)
plot(distr_cc, col=mypalette(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(europe, add=T)
mtext("Climate scenario RCP85 HD ", 3,-2.2)
mtext("Cicendia filiformis ", 3,-1.2, font=3)
mtext("WGS84 ", 1,-1.2, adj=1)
dev.off()

### Q: Is this expected change of distribution reasonably from an ecological perspective?
### Q: Schould we add a legend? If yes, what's the problem and how does a legend make sense?

### Uncertainty: Derive from validation, e.g. REPEATED DATA SPLITTING!!!

