# --------------------------------- Figure 4 --------------------------------- #

library(DataFission)
library(datathin)
library(ggplot2)
library(pbapply)
library(dplyr)
library(latex2exp)
library(patchwork)
theme_set(theme_bw())

#-- Functions 
source("utils.R")

kernel <- function(u, h) {
  a <- h *sqrt(3)
  return(ifelse(abs(u) < a, 0.5/a, 0))
}

weightvar <- function(x,obs, h = NULL){
  if (is.null(h)) {
    h <- sd(x) * (4/(3*n))^(1/5)
  }
  weigth <- kernel(x-obs, h)
  return(modi::weighted.var(obs, weigth))
}

local_var <- function(x){
  n <- length(x)
  h0 <- ks::hpi(x, deriv.order = 1)
  pilot <- ks::kde(x, h = h0, eval.points=x)$estimate
  lambda <- exp(mean(log(pilot)))
  hi <- h0*sqrt(lambda/pilot)
  variance_wk <- sapply(1:n, FUN = function(i){weightvar(x[i], obs = x, h = hi[i])})
  return(variance_wk)
}


sim_fun_power <- function(seed, delta=20, n=100, sd = 1, tau = .4){
  set.seed(250124*seed)
  X <- c(rnorm(n/2, mean = 0, sd = sd),
         rnorm(n/2, mean = delta, sd = sd))
  variance <- local_var(X)
  Z <- sapply(variance, function(s){rnorm(1, mean = 0, sd = sqrt(s))})
  fX  <- X + tau*Z 
  gX <- X - (1/tau)*Z 
  
  cl <- km_fun(fX, K = 2)
  df <- data.frame(gX = gX, 
                   Cluster = as.factor(cl))
  
  res.ttest <- t.test(gX~cl)$p.value
  fiss_glob <- data_fission(as.matrix(X), tau = tau)
  cl_glob <- km_fun(fiss_glob$fX, K=2)
  res.ttestglob <- t.test(fiss_glob$gX~cl_glob)$p.value
  return(data.frame(pval = c(res.ttest, res.ttestglob),
                    ARI = c(mclust::adjustedRandIndex(cl, rep(1:2, each = n/2)),
                            mclust::adjustedRandIndex(cl_glob, rep(1:2, each = n/2))),
                    Variance = c("Local", "Global")))
}

apply_sim_delta <- function(seed, delta_grid, n=100, sd = 1, tau = .4){
  res_delta <- lapply(delta_grid, function(d){
    res_temp <- sim_fun_power(seed = seed, delta = d, sd = sd, tau = tau)
    res_temp$delta <- d
    return(res_temp)})
  temp <- do.call("rbind.data.frame", res_delta)
  return(temp)
}


apply_sim_delta_tau <- function(seed, delta_grid, n=100, sd = 1, tau_grid){
  res_tau <- lapply(tau_grid, function(t){
    res_temp <- apply_sim_delta(seed = seed, delta_grid = delta_grid, n = n, sd = 1, tau = t)
    res_temp$tau <- t
    return(res_temp)
  })
  temp <- do.call("rbind.data.frame", res_tau)
  return(temp)
}

n <- 100
delta_grid <- c(3.5, 5, 8, 10, 20)
tau <- seq(0.1, 3, length.out = 50)
sigma <- 0.5
nsimu <- 100

res_sim <- pblapply(1:nsimu, function(ns){
  res <- apply_sim_delta_tau(seed = ns, delta_grid = delta_grid, n = n, sd = sigma, tau_grid = tau)
  return(res)
}, cl = 5)

df <- do.call(rbind.data.frame,res_sim) %>%
  group_by(delta, tau, Variance) %>%
  summarise(Power = mean(pval < 0.05), 
            MeanARI= mean(ARI),
            sdARI = sd(ARI))

plot_typeI <- df %>% 
  ggplot() + 
  aes(x=tau, y = Power, colour = as.factor(delta), linetype = Variance) +
  geom_line(size = 1.4) +
  NULL

plot_ARI <- df %>% 
  ggplot() + 
  aes(x=tau, y = MeanARI, colour = as.factor(delta), linetype = Variance) +  # geom_point(size = 3) +
  geom_line(size = 1.4) +
  NULL

  