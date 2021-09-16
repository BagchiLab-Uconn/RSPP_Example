# initialize ####

rm(list = ls())

require(RSPPlme4)
require(dplyr)

# load dataset ####

data <- spatstat.data::lansing

# these data start as a ppp (planar point pattern), which is an easier starting point than most datasets
# so, first, this code will convert that ppp to a dataframe, which is a more standard starting point
df_data <- as.data.frame(data)

# as "misc" is not a single tree species, we will not use it in this example
df_data <- df_data[!df_data$marks %in% "misc",]

# organize dataset ####

# in this example, we will analyze two things:
# (1) how different is the spatial clustering of trees between species?
# (2) how can using replicated spatial point patterns reduce overall sampling effort?

# first, since this lansing plot is a single, large plot, 
# we will randomly sample 20 smaller plots within, comprising ~20% of the total site

sampled_SPs <- 20 # number of subplots to sample
percentArea <- .4 # approx percent of total area to sample
totalPlots <- sqrt(sampled_SPs/percentArea) %/% 1 # is actually the sqrt of total plots

# these will the x and y boundaries of the sub-plots
xy <- seq(0, 1, length.out = totalPlots + 1)

# these are the specific boundaries of the specific subPlots
subPlots <- data.frame(subPlotID = (1:(totalPlots)^2), 
                    xmin = rep(xy[-totalPlots], totalPlots),
                    ymin = rep(xy[-totalPlots], each = totalPlots))

# this determines the subPlot of each tree 

df_data$xmin <- df_data$x %/% xy[2] * xy[2]
df_data$ymin <- df_data$y %/% xy[2] * xy[2]
divided_data <- merge(df_data, subPlots, by = c("xmin","ymin"))
divided_data$x <- divided_data$x - divided_data$xmin
divided_data$y   <- divided_data$y - divided_data$ymin
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
  pppx <- ppp(x=dat$x, y = dat$y, window = owin(x=c(0,xy[2]), y=c(0, xy[2])))
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

# stall out here due to spatstat issue







