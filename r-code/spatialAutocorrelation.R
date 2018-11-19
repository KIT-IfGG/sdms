# packages used for the data generation
library(raster)
library(colorRamps) # for some crispy colors
library(vegan) # will be used for PCNM
library(ncf) # correlogramm (there are alternatives!)
library(spdep) # neighbours

### SIMULATE DATA ####
# empty matrix and spatial coordinates of its cells
side <- 30
my.mat <- matrix(NA, nrow=side, ncol=side)
x.coord <- rep(1:side, each=side)
y.coord <- rep(1:side, times=side)
xy <- data.frame(x.coord, y.coord)

### Environmental effect ####
# Create environemntal variable
#env.value <- rnorm(side*side, 0, 2)
env.value <- 0.5 * xy$x.coord - 0.007 * xy$y.coord + 0.002 * xy$y.coord^2 

# Write value into matrix
env.mat <- my.mat
env.mat[] <- env.value
r.env <- raster(env.mat)

# Calculate response variable based on env
coeffs <- c(3, 0.9, 0.2)
z.value_lin <- coeffs[1] + coeffs[2] * env.value + coeffs[3] * env.value^2 + rnorm(side*side, 0, 1)
my.mat[] <- z.value_lin
r.response <- raster(my.mat)

### Plot environemental variable
x11(width=8, height=8)
#pdf("env.pdf", width=10, height=5)
par(mfrow=c(2,2))
plot(r.env, axes=F, col=matlab.like(20), main="Environmental variable")

### Plot response variable ####
plot(r.response, axes=F, col=matlab.like(20), main="Response variable")

### Linear regression ####
dat_lin <- data.frame(env.value, z.value_lin)
m <- lm(z.value_lin ~ env.value + I(env.value^2), dat_lin)
round(coef(m), 1)
coeffs
### -> Parameters used to create the data can be found by lm.

### Spatial effect ####
# all paiwise euclidean distances between the cells
xy.dist <- dist(xy)
# PCNM axes of the dist. matrix (from 'vegan' package)
pcnm.axes <- pcnm(xy.dist)$vectors
# using 8th PCNM axis as my atificial z variable
z.space <- pcnm.axes[,8]*100 + rnorm(side*side, 0, 1)
z.space <- scale(z.space, center=TRUE, scale=FALSE)
# plotting the artificial spatial effect
my.mat[] <- z.space
r.space <- raster(my.mat)
plot(r.space, axes=F, col=matlab.like(20))

### Neighbours and their correlation ####
xy.neigh <- dnearneigh(as.matrix(xy), d1=0, d2=4, longlat = FALSE)
str(xy.neigh, 2)

# 1st order neighbour of element 1
xy.neigh[[1]][1]
# 1st order neighbour of element 3
xy.neigh[[3]][1]
# 3rd order neighbour of element 1
xy.neigh[[1]][3]

neigh_1st <- unlist(lapply(xy.neigh, function (x) x[1]))
cor(my.mat[1:length(my.mat)], my.mat[neigh_1st], method="pearson")
cor(my.mat[1:length(my.mat)], my.mat[neigh_1st], method="spearman")
cor(my.mat[1:length(my.mat)], my.mat[neigh_1st], method="kendall")
## -> Neighbours are corelated

# Plot correlation
spearman.cor <- numeric()
for (i in 1:10) spearman.cor[i] <- cor(my.mat[1:length(my.mat)], my.mat[unlist(lapply(xy.neigh, function (x) x[i]))], method="spearman")
plot(1:length(spearman.cor), spearman.cor, xlab="Order of Neighbourhood", ylab="Pearson r", type="h")

### Correlogram ####
ncf.cor <- correlog(x.coord, y.coord, z.space, increment=2, resamp=10)
plot(ncf.cor)

### Create data with environmental and spatial effects ####
z.value <- z.space + z.value_lin
my.mat[] <- z.value
r.z <- raster(my.mat)

#pdf("z.pdf", width=5, height=5)
x11()
plot(r.z, axes=F, col=matlab.like(20))
dev.off()

### MODELS ####
dat <- data.frame(z.value = z.value, env.value = env.value)

my.lm <- lm(z.value ~ env.value + I(env.value^2), data=dat)
summary(my.lm)

x11()
plot(residuals(my.lm) ~ env.value)
abline(0,0, col="red")

### Plot residuals in space
resid.mat <- my.mat
resid.mat[] <- residuals(my.lm)
r <- raster(resid.mat)
plot(r, axes=F, col=matlab.like(20))
### - > Residuals show clear spatial autocorrelation!

### Correlogram of residuals
cor.resids <- correlog(x.coord, y.coord, residuals(my.lm), increment=2, resamp=10)
plot(cor.resids)

### Better model
my.lm.space <- lm(z.value ~ env.value + I(env.value^2), data=dat, offset=z.space)   ### Offset is not really a way to deal with autocorrelation! better: coorelation structure in a marginal model.
summary(my.lm.space)

### Compare coefs with real values
round(coef(my.lm),2)
round(coef(my.lm.space),2)
coeffs

cor.resids.space <- correlog(x.coord, y.coord, residuals(my.lm.space), increment=2, resamp=10)

x11()
plot(cor.resids.space)

resid.mat[] <- residuals(my.lm.space)
r <- raster(resid.mat)
plot(r, axes=F, col=matlab.like(20))
### -> No spatial autocorrelation anymore :-).

### Alternative packages for spatial autocor ####
# library(pgirmess)
# pgi.cor <- correlog(coords=xy, z=z.value, method="Moran", nbclass=21)
# 'nb' - neighbourhood of each cell
# r.nb <- dnearneigh(as.matrix(xy), d1=0.5, d2=1.5)
# 'nb' - an alternative way to specify the neighbourhood
# r.nb <- cell2nb(nrow=side, ncol=side, type="queen")
# sp.cor <- sp.correlogram(r.nb, z.value, order=15, method="I", randomisation=FALSE)

### Marginal models ####
### Use correaltion structure for considering spatial autocorrelation of residuals.
library(lme4)
library(nlme)
dat <- cbind.data.frame(dat, xy)
my_gls <- gls(z.value ~ env.value + I(env.value^2), data=dat, correlation=corGaus(form= ~ x.coord + y.coord))   ### Offset is not really a way to deal with autocorrelation! Better: correlation structure

summary(my_gls)
### Compare coefs with real values
round(coef(my_gls),1)
coeffs

### Correlagramm
cor.resids <- correlog(x.coord, y.coord, residuals(my_gls), increment=2, resamp=10)
plot(cor.resids)

### Plot residuals in space
resid.mat <- my.mat
resid.mat[] <- residuals(my_gls)
r <- raster(resid.mat)
plot(r, axes=F, col=matlab.like(20))

### This correlation is to strong/compex for a simply correlaiton structure. There are other methods more suitable for this case.

### Spatial Eigenvectors (SEVM) ####
# Correcting for spatial autocorrelation using Spatial Eigenvectors (SEVM)
# Various other methods explained in the Dormann et al. 2007 paper

# Identify neighbors
xy.neigh <- dnearneigh(as.matrix(xy), d1=0, d2=10, longlat = FALSE)
# Checken ob irgendein Punkt nicht angebunden ist
summary(xy.neigh)
xy.neigh_w <- nb2listw(xy.neigh)
res <- ME(z.value ~ env.value, dat, family = gaussian, listw = xy.neigh_w)
dat <- cbind.data.frame(dat, fitted(res))
m <- lm(z.value ~ env.value + I(env.value^2) + vec4 + vec5 + vec10 + vec7, data = dat)
m
coeffs
