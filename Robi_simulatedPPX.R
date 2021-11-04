#devtools::install_github("Bagchilab-Uconn/RSPPlme4")
library(RSPPlme4)
library(spatstat)
library(ggplot2)
library(tidyverse)

# Data simulation
set.seed(1234)
n <-  100

dat <- hyperframe(x = runif(n = n, min = 0, max = 10),
                  f = sample(c("a", "b"), size = n, replace = TRUE),
                  gr = rep(1:10, each = n/10))

lin_pred1 <-  model.matrix(~ f + x, data = as.data.frame(dat, warn = FALSE)) %*%
  c(0.05, 0.05, 0.01)

dat$ppx1 <- lapply(lin_pred1, function(sigma)
                   rThomas(kappa =runif(n = 1, min = 3, max = 10),
                           mu = 5, sigma = sigma))

ranef <- exp(rnorm(n = 10, mean = 0, sd = 0.2) )

lin_pred2  <- lin_pred1 * ranef[dat$gr]

dat$ppx2 <- lapply(lin_pred2, function(sigma)
  rThomas(kappa =runif(n = 1, min = 3, max = 10),
          mu = 10, sigma = sigma))

# Models
mod1.1 <- klm(ppx1 ~ 1, hyper = dat, r = seq(0, 0.1, 0.02), 
            correction = "border", weights_type = "nx_A")

mod1.1_ci <- confint(mod1.1, nboot=500, level = 0.95, iseed = 4321)
print(mod1.1_ci)

plot(mod1.1_ci) + ggthemes::theme_tufte() + #scale_y_continuous(trans = "sqrt") +
  geom_function(fun = function(x) pi * x^2, linetype = "dotted") + 
  labs(x = "Distance (r)", y = "K(r)")

# add covariates
mod1.2 <- klm(ppx1 ~ f + x, hyper = dat, r = seq(0, 0.1, 0.01),
            correction = "border", weights_type = "nx_A")

mod1.2_ci <- confint(mod1.2, nboot = 500, level = 0.99, iseed = 9876) 

plot(mod1.2_ci) + ggthemes::theme_tufte() 


nd1 <- expand.grid(x = c(1, 5, 10),  f = c("a", "b"))
mod1.2_ci <- confint(mod1.2, nboot = 500, level = 0.99, iseed = 9876, 
                   newdata=nd1) 

preds1.2  <- as.data.frame.table(mod1.2_ci$predictions)
preds1.2[, c("x", "f")] <- nd1[as.numeric(preds1.2$Var1), c("x", "f")]
preds1.2 <- pivot_wider(preds1.2, names_from = Var3, values_from = Freq) %>%
  rename("distance" = "Var2") %>%
  select(-Var1) %>% mutate(distance = as.numeric(distance))

ggplot(preds1.2, aes(x = distance, y = est, ymin = lwr, ymax = upr, 
                     group = as.factor(x))) +
  geom_ribbon(alpha = 0.2, aes(fill = as.factor(x))) + 
  geom_path(aes(color = as.factor(x))) +
  facet_wrap(~f) + 
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Set2") +
  ggthemes::theme_tufte() + scale_y_continuous(trans = "L") + 
  labs(x = "distance (r)", y = expression(italic(L(r))), color = "x", fill = "x")


mod2 <- klmer(ppx2 ~ 1 + f + x +  (1| gr), hyper = dat, r = seq(0, 0.1, 0.01),
              correction = "border", weights_type = "nx_A")

mod2_ci <- confint(mod2, nboot = 500, ncore = 16)

print(mod2_ci)  
plot(mod2_ci) + ggthemes::theme_tufte() +
  labs(x = "Distance (r)", y = "K(r)")

mod2_ci <- confint(mod2, nboot = 99, ncore = 16, 
                   newdata  = nd1)

preds2  <- as.data.frame.table(mod2_ci$predictions)
preds2[, c("x", "f")] <- nd1[as.numeric(preds2$Var2), c("x", "f")]
preds2 <- pivot_wider(preds2, names_from = Var3, values_from = Freq) %>%
  rename("distance" = "Var1") %>%
  mutate(distance = as.numeric(as.character(distance)))

ggplot(preds2, aes(x = distance, y = est, ymin = lwr, ymax = upr, 
                   colour = as.factor(x), fill = as.factor(x))) + facet_wrap(~f) + 
  geom_ribbon(alpha = 0.1, color = NA, aes(group = as.factor(x))) +geom_path() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Set2") +
  scale_y_continuous(trans = "L") + 
  labs(x = "distance (r)", y = expression(italic(L(r))), color = "x", fill = "x") +
  ggthemes::theme_tufte() 

