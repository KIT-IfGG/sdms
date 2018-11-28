library(gjam)
library(vegan)
library(classInt)

### Das ist ein Beispiel  ### das ist nicht richtig

### Read data ####
veg <- read.table("data/kreuz.txt")
env <- read.table("data/env.txt")

### Check data ####
summary(veg)
summary(env)
table(rownames(veg) == rownames(env))
# -> All ok!

### Data description / data mining ####

### Data transformation
veg_trans <- decostand(veg, "max")

### Ordination
ordi <- metaMDS(veg_trans, distance = "bray", k = 2, trace = TRUE)   ### NMDS
ef <- envfit(ordi, env[,c("Jdatum", "DDSUMLuft", "WG.min10")])
oeh <- ordihull(ordi, env[,"Lokalitaet"])


### Plot ordination results
ordiplot(ordi, choices=c(1,2), type = "n")
orditorp(ordi, choices=c(1,2), display="species", cex=1, priority=colSums(veg_trans), col = "black")
points(ordi, choices=c(1,2), display="sites", col = rainbow(nlevels(env$Plot), alpha = 0.5)[as.numeric(env$Plot)], cex=1, pch=16)
plot(ef, add=T, col = "darkblue")
ordihull(ordi, env[,"Lokalitaet"], label = TRUE, col = "darkblue", cex = 1.5)

### JSDM ####
dat <- list()
dat$formula <- as.formula("~ DDSUMLuft + WG.min10")
dat$xdata <- env[,c("DDSUMLuft", "WG.min10")]
dat$ydata <- veg

ml  <- list(ng = 1000, burnin = 1000, typeNames = "CA")

out <- gjam(dat$formula, dat$xdata, dat$ydata, modelList = ml)
saveRDS(out, "data/jsdm_out.rds")

### Mini - assessment ####
species_fit <- gjamPredict(output = out, newdata = list(xdata = dat$xdata))
corres <- numeric()
for(i in 1:ncol(veg)) {
  sp <- colnames(veg)[i]
  corres[i] <- round(cor(veg[,sp], species_fit$sdList$yMu[,sp], method = "spearman"), 3)
}

hist(corres)
mean(corres)

### Response curves ####
n <- 20
newdata <- expand.grid(DDSUMLuft = seq(min(dat$xdata$DDSUMLuft), max(dat$xdata$DDSUMLuft), len = n), WG.min10 = seq(min(dat$xdata$WG.min10), max(dat$xdata$WG.min10), len = n))

species_pred <- gjamPredict(output = out, newdata = list(xdata = newdata))
#saveRDS(species_pred, file="results/species_preds.rds")
newdata <- cbind.data.frame(newdata, species_pred$sdList$yMu)

### Select species
sp <- "Cicefili"

ci <- classIntervals(newdata[,sp], n)
spcols <- findColours(ci, colorRampPalette(c("grey", "green"))(n))
xdiff <-diff(unique(newdata$DDSUMLuft))[1]
ydiff <-diff(unique(newdata$WG.min10))[1]


symbols(x = newdata$DDSUMLuft, y = newdata$WG.min10, rectangles = matrix(rep(c(xdiff, ydiff), nrow(newdata)), ncol=2, byrow=T), bg = spcols, fg = spcols, inches = F, xlab = "DD Sum Luft", ylab = "min Wassergehalt 10 Tage")
points(x = dat$xdata$DDSUMLuft, y = dat$xdata$WG.min10, col=ifelse(veg[,sp] > 0, "black", "white"), pch=16, cex = 1.4)
mtext(sp, 3, -1.5)

### Loop for all species ####
pdf("figures/jsdm_species_response.pdf", height = 9, width = 6)
par(mfrow=c(3,2))
for(i in 3:ncol(newdata)) {
  sp <- colnames(newdata)[i]  
  ci <- classIntervals(newdata[,sp], n)
  spcols <- findColours(ci, colorRampPalette(c("grey", "darkgreen"))(n))
  xdiff <-diff(unique(newdata$DDSUMLuft))[1]
  ydiff <-diff(unique(newdata$WG.min10))[1]
  symbols(x = newdata$DDSUMLuft, y = newdata$WG.min10, rectangles = matrix(rep(c(xdiff, ydiff), nrow(newdata)), ncol=2, byrow=T), bg = spcols, fg = spcols, inches = F, xlab = "DD Sum Luft", ylab = "min Wassergehalt 10 Tage")
  points(x = dat$xdata$DDSUMLuft, y = dat$xdata$WG.min10, col=ifelse(veg[,sp] > 0, "black", "lightgrey"), pch=16, cex = 1.2)
  mtext(paste0("Spearman r = ", corres[i], " "), 1, -1.5, adj = 1)
  mtext(sp, 3, -1.5)
}
dev.off()  


