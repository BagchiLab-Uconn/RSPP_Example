# initialize ####

rm(list = ls())

require(spatstat)
require(RSPPlme4)
require(tidyverse)

# load dataset 
data <- spatstat.data::lansing

# as "misc" is not a single tree species, we will not use it in this example
data <- subset(data, marks != "misc", drop = T)

# these data start as a ppp (planar point pattern), which is an easier starting point than most datasets
# so, first, this code will convert that ppp to a dataframe, which is a more standard starting point
df_data <- as.data.frame(data)

# organize dataset ####

# in this example, we will analyze two things:
# (1) how different is the spatial clustering of trees between species?
# (2) how can using replicated spatial point patterns reduce overall sampling effort?

# first, since this lansing plot is a single, large plot, 
# we will randomly sample 20 smaller plots within, comprising ~20% of the total site

sampled_SPs <- 20 # number of subplots to sample
percentArea <- .6 # approx percent of total area to sample
totalPlots <- sqrt(sampled_SPs/percentArea) %/% 1 # is actually the sqrt of total plots

# these will the x and y boundaries of the sub-plots
xy <- seq(0, 1, length.out = totalPlots + 1)

# these are the specific boundaries of the specific subPlots
subPlots <- data.frame(subPlotID = (1:(totalPlots)^2), 
                    xmin = rep(xy[-length(xy)], totalPlots),
                    ymin = rep(xy[-length(xy)], each = totalPlots))

# this determines the subPlot of each tree 

df_data$xmin <- df_data$x %/% xy[2] * xy[2]
df_data$ymin <- df_data$y %/% xy[2] * xy[2]
divided_data <- merge(df_data, subPlots, by = c("xmin","ymin"))
divided_data$x <- (divided_data$x - divided_data$xmin) * 100 / xy[2]
divided_data$y   <- (divided_data$y - divided_data$ymin) * 100 / xy[2]
divided_data$xmin <- NULL
divided_data$ymin <- NULL

# sample 20 plots

set.seed(11)
ourPlots <- (divided_data[divided_data$subPlotID %in% sample(1:max(divided_data$subPlotID), sampled_SPs),])


# convert the data to a hyperframe
grb <- ourPlots %>% 
  group_by(subPlotID, marks) %>%
  summarise()
hf <- as.hyperframe(grb)

# convert each plant-subplot into its own ppp

getPPP <- function(dat){
  dat <- dat[!duplicated(dat[, c("x", "y")]),]
  pppx <- ppp(x=dat$x, y = dat$y, window = owin(x=c(0,100), y=c(0, 100)))
}
spp <- lapply(split(ourPlots, 
                    f = list(ourPlots$subPlotID, ourPlots$marks), 
                    drop=TRUE),
              getPPP)

hf$pppx <- spp

# remove any plant-subplot ppp that has 3 or fewer points

dat <- hf[sapply(hf$pppx, npoints) > 3,]

# create K estimate 
r <-  0:25
dat$K <-  lapply(dat$pppx, Kest, r = r, correction = "border")
dat$wts <- lapply(dat$pppx, kfuncWeightsCalc, r=r,
                  correction="border", type="nx_A")

# check none of the weights are 0
which(sapply(dat$wts, function(x) any(x ==0)))

# remove those rows
dat <- dat[!(sapply(dat$wts, function(x) any(x ==0))),]

# run model
mod_sm <- klmerHyper(formula = K ~ 1  + marks + (1|subPlotID),
                     r=r, hyper=dat,  correction='border', weights = wts,
                     na.action="na.omit")

# plot models

plot(mod_sm)
a <- anova(mod_sm, dists = 1:25, term = "marks", nboot = 10)
a

preddat <- expand.grid(marks = factor(c("blackoak", "hickory", "maple", "redoak", "whiteoak")))
Xmat <- model.matrix(~marks, data=preddat)  
sapply(mod_sm, function(x) VarCorr(x))
summary(mod_sm$`1`)
mod_sm_cis <- confint(mod_sm, level=0.95, lin_comb=Xmat, nboot=10, ncore=4, iseed=1234)
plot(mod_sm_cis)


klmerci2plot(mod_sm_cis, preddat) %>%
  ggplot(aes(x=distance, y=est, fill= marks,
             group= marks,
             ymin = lcl_pred, ymax = ucl_pred)) +
  geom_ribbon(alpha=0.3) + geom_path(aes(colour= marks)) +
  facet_wrap(aes(marks))

# full plot analysis
dd <- divided_data %>%
  split(f = divided_data$marks, 
        drop = T)
fp <- lapply(dd, function(dat) ppp(x=dat$x, y = dat$y, window = owin(x=c(0,100), y=c(0, 100))))

ks <- lapply(fp, Kest, correction = "border") 
par(mfcol = c(3,2))
lapply(ks, plot)
ests <- lapply(fp, envelope)
lapply(ests, plot)


# now what if 2 species are thinned, and two species are thickened?
dev.off()
plot(data)





byPlant <- split(data, f = data$marks, drop = T)


kdens <- density(byPlant, at = "points", sigma = 0.05)
head(kdens)
hist(kdens[[3]])
maple <- byPlant[[3]]
maple$marks <- factor(maple$marks)

plot(maple, cols = kdens[[3]])

library(ggplot2)

ggplot(mapping = aes(x = maple$x, y = maple$y)) +
  geom_point(aes(color = kdens[[3]]), size = 2) +
  scale_color_viridis_c() +
  theme_bw()

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
i <- 3
hist(range01(log(density(byPlant[[i]])/max(density(byPlant[[i]])))))



# thinned <- lapply(3:4, function(i) rthin(byPlant[[i]], P = 1 - plogis(qlogis(0.25) + 10 * range01(kdens[[i]])^2)))
# par(mfcol = c(2,2))
# estnk <- lapply(thinned, envelope, correction = "border")
# lapply(estnk, plot)
# lapply(thinned, plot)
# lapply(thinned, plot)
# lapply(byPlant[3:4], plot)

dev.off()
par(mfcol = c(2,2))

lapply(1:5, function(i){
  maple <- byPlant[[i]]
  mm <- as.matrix(dist(coords(maple)))
  mm[lower.tri(mm)] <- 0
  
  thold <- .04
  keepers <- apply(mm, 1, function(x, told = thold){
    y <- x[x != 0]
    if(length(y) != 0){
      if(min(y) < told)
        keep <- 0
      else
        keep <- 1
    }
    else
      keep <- 1
    return(keep)
  })
  
  tmap <- rthin(maple, P = keepers)
  
  
  plot(envelope(maple, correction = "border"))
  plot(envelope(tmap, correction = "border"))
  
  plot(maple)
  plot(tmap)
})


redoak <- byPlant[[4]]
newh <- rpoint(redoak$n, f = range01(density(redoak)^2))
newh$marks <- as.factor(rep("redoak", newh$n))
newh$markformat <- "vector"
hs <- superimpose(redoak, newh)

plot(envelope(hs, correction = "border"))
plot(envelope(redoak, correction = "border"))

plot(hs)
plot(redoak)




clust <- function(x,y){
  xx <- x < .25 & y < .25
  return(as.numeric(xx))
}

even <- data.frame(x = rep(seq(0,1, by = .1), 11), y = rep(seq(0,1, by = .1), each =11))
p <- ppp(even$x, even$y, window = owin(c(0,1), c(0,1)))
pj <- rjitter(p, radius = .01)
plot(pj)
plot(envelope(pj, correction = "border"))

plot(rpoint(100, f = clust))


rr <-rpoint(100, f = clust)
plot(envelope(rr, correction = "border"))

uu <- runifpoint(100)
plot(envelope(uu, correction = "border"))

