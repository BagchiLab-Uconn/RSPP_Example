# initialize ####

rm(list = ls())

library(spatstat)
library(sp)
library(tidyverse)
library(ggthemes)
library(patchwork)


set.seed(123)
# make a uniform ppp ####

# what is the sqaure root of the number of points you want?
sqrtNumPoints <- 10
win <- owin(c(0, 100), c(0, 100))
even <- data.frame(x = rep(seq(1/sqrtNumPoints/2 , 100 - 1/sqrtNumPoints/2, 
                               length.out = sqrtNumPoints), sqrtNumPoints),
                   y = rep(seq(1/sqrtNumPoints/2 , 100 - 1/sqrtNumPoints/2, 
                               length.out = sqrtNumPoints), each =sqrtNumPoints))

p <- ppp(even$x, even$y, window = win)
uniform <- rjitter(p, radius = 3)
plot(uniform)

# make a clustered ppp ####
# make ppp 
clustered <- rThomas(sqrtNumPoints/1e4, scale = 3, mu =  10, 
                           win=owin(c(0, 100), c(0, 100)))
plot(clustered)

# make a random ppp ####
rando <- runifpoint(sqrtNumPoints^2, win = owin(c(0,100), c(0,100)))

# make k estimates
r <- 0:25
methods(class=class(kr)[1])
ku <- Kest(uniform, r = r, correction = "border")
kr <- Kest(rando, r = r, correction = "border")
kc <- Kest(clustered, r = r, correction = "border")

nullmod <- envelope(rando, correction = "border", nsim=1e3)
nullmod <- data.frame(nullmod)

# make k estimate dataframe
curves <- data.frame(
  type = factor(rep(c("inhibited", "random", "clustered"), each = length(r)),
                levels = c("clustered", "random", "inhibited" )),
  r = rep(r, 3),
  Kr = c(ku$border, kr$border, kc$border))

## Make palette
pal <- RColorBrewer::brewer.pal(n=3, name ="Dark2")
K2L <- function(x) sqrt(x/pi)
# make k estimate plot ####
p4 <- ggplot(data = curves, aes(x = r)) +
  xlim(0,25) +
  geom_function(fun = function(x) pi * x^2, color = "black", linetype = "dashed") +
  geom_ribbon(data = nullmod, aes(ymin = lo, ymax = hi), alpha = .15) +
  geom_line(aes(color = type, y = Kr), size = 1) + 
  theme_tufte() +
  scale_colour_brewer(palette="Dark2") +
  scale_y_continuous(breaks=seq(0, 2000, 500)) +
  labs(y ="K(r)", colour = NULL ) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

# make L estimate plot
p5 <- ggplot(data = curves, aes(x = r)) +
  xlim(0,25) +
  geom_function(fun = function(x)  x, color = "black", linetype = "dashed") +
  geom_ribbon(data = nullmod, aes(ymin = K2L(lo), ymax = K2L(hi)), alpha = .15) +
  geom_line(aes(color = type, y = K2L(Kr)), size = 1) + 
  theme_tufte() +
  scale_colour_brewer(palette="Dark2") +
  labs(y = "L(r)", x = "distance, r", colour = NULL )
  

# make ppp plots ####
plotit <- function(p, clr){ ggplot(data = coords(p), aes(x = x, y = y)) +
    geom_point() +
    theme_bw() +
    coord_equal() +    
    theme(legend.position = "none",
          panel.border = element_rect(color = clr, fill = NA),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
}

p1 <- plotit(clustered, pal[1] )
p2 <- plotit(rando, pal[2])
p3 <- plotit(uniform, pal[3])


# make full panel patchwork ####
p1 <- p1 + labs(tag="A")
p4 <- p4 + labs(tag ="B")
p5 <- p5 + labs(tag ="C")


((p1 / p2 / p3) | (p4/p5)) + 
  plot_layout(widths = c(1.5, 2), guides = "collect") &
    theme(legend.position = "right")




