library(spdep)
library(akima)
library(rgl)
library(ncf)
library(nlme)

### Simulate data ####
### Space ####
N <- 200
x.coord <- runif(N,0,100)
y.coord <- runif(N,0,100)
points <- cbind(x.coord,y.coord)
plot(y.coord ~ x.coord, data=points, cex = 0.6, col = "blue", pch = 16)

### Space to environment ####
coeffs.env <- c(0, 0.06, 0.03)
points <- as.data.frame(points)
env <- coeffs.env[1] + coeffs.env[2] * points$x.coord + coeffs.env[3] * points$y.coord

### Environemnt to response variable ####
coeffs <- c(alpha = 2, beta = 0.8, range = 3.5)   ### Regression parameters to be estimated!
z.env <-  rnorm(N, coeffs[1] + coeffs[2] * env, coeffs[2]/10)

### Random variables spatially autocorrelated ####
z.space <- cbind(rmvn.spa(x=points$x.coord, y=points$y.coord, p=coeffs[3], method="gaus")) * sd(z.env)
hist(z.space); mean(z.space); range(z.space)

### Just to show the spatial pattern! It is interpolated!
surv.space <- interp(points$x.coord, points$y.coord, z.space[,1], linear=TRUE)
image(surv.space$x, surv.space$y, surv.space$z, main = "Autocorrelated random variable")
contour(surv.space$x, surv.space$y, surv.space$z, add=T)
mean(z.space)

### Add spatial autocorr "error" to deterministic response to env ####
z <- z.env + z.space

points <- cbind.data.frame(x.coord, y.coord, z, env)

### Plot, interpolated!
surv <- interp(points$x.coord, points$y.coord, points$z, linear=TRUE)
image(surv$x, surv$y, surv$z, main = "Response variable")
contour(surv$x, surv$y, surv$z, n = 5, add = T)
points(points$x, points$y, cex = 0.6, col = "blue", pch = 16)

### Model ####
#m1 <- lm(z ~ env, data = points)
m1 <- gls(z ~ env, data = points)
m2 <- gls(z ~ env, data = points, correlation = corGaus(form = ~ x.coord + y.coord))

### Compare fitted parameters
rbind(m1 = round(coef(m1),2), m2 = round(coef(m2),2), coeffs = coeffs[1:2])

BIC(m1)
BIC(m2)

summary(m1)
summary(m2)

points$residuals <- resid(m2)

surv.resids <- interp(points$x.coord, points$y.coord, points$residuals, linear=TRUE)
image(surv.resids$x, surv.resids$y, surv.resids$z, main = "Residuals")

ncf.cor <- correlog(x = points$x.coord, y = points$y.coord, z = points$residuals, increment = 2, resamp = 100)
#ncf.cor$n
plot(ncf.cor$mean.of.class, ncf.cor$correlation, type = "b", pch = 16, col=ifelse(ncf.cor$p < 0.05, "red", "black"), xlab = "Distance", ylab = "Spatial autocorrelation", xlim=c(0, 40))
abline(h = 0, col = "blue")
abline(v = coeffs[3], col = "black")
text(coeffs[3]-1, max(ncf.cor$correlation)/2, "Range of autocorrelation", srt = 90, lty=2)
legend("bottomright", legend=c("Sign.", "non sign."), col=c("red", "black"), pch=16)

