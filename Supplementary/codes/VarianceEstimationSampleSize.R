#------------------------- Supplementary materials ----------------------------#
#------------- Impact of the sample size on variance estiamtion ---------------#

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


sim_fun <- function(seed, delta=20, n, sd = 1, tau = 0.4){
  set.seed(250124*seed)
  X <- c(rnorm(n/2, mean = 0, sd = sd),
         rnorm(n/2, mean = delta, sd = sd))
  variance <- local_var(X)
  return(data.frame(Variance_hat = mean(variance)))
}

apply_sim_delta <- function(seed, delta_grid, n, sd = 1, tau = 0.4){
  res_delta <- lapply(delta_grid, function(d){
    res_temp <- sim_fun(seed = seed, delta = d, sd = sd, n=n)
    return(res_temp)})
  temp <- do.call("rbind.data.frame", res_delta)
  temp$delta <- delta_grid
  return(temp)
}

apply_sim_delta_sigma <- function(seed, delta_grid, n, sd_grid, tau = 0.4){
  res_sigma <- lapply(sd_grid, function(s){
    res_temp <- apply_sim_delta(seed = seed, delta_grid = delta_grid, n = n, sd = s, tau = .4)
    res_temp$sigma <- s
    return(res_temp)
  })
  temp <- do.call("rbind.data.frame", res_sigma)
  return(temp)
}


apply_sim_delta_sigma_n <- function(seed, delta_grid, n_grid, sd_grid, tau){
  res_n <- lapply(n_grid, function(samp_size){
    res_temp <- apply_sim_delta_sigma(seed = seed, delta_grid=delta_grid, n=samp_size, sd_grid = sd_grid, tau = tau)
    res_temp$SampSize <- samp_size
    return(res_temp)
  })
  temp <- do.call("rbind.data.frame", res_n)
  return(temp)
}

#--Param 
n <- c(5, 10, 25, 50, 100, 250, 500)
delta_grid <- c(seq(0, 3, length.out = 25), seq(3.5, 100, length.out = 25))
tau <- .4
sigma <- c(0.1, 0.5, 1, 2) 
nsimu <- 100

res_sim <- pblapply(1:nsimu, function(ns){
  res <- apply_sim_delta_sigma_n(seed = ns, delta_grid = delta_grid, n_grid = n, sd_grid = sigma, tau = tau)
  return(res)
}, cl = 5)

df <- do.call(rbind.data.frame,res_sim) %>%
  group_by(delta, sigma, SampSize) %>%
  summarise(Variance= mean(Variance_hat),
            sdVariance_hat = sd(Variance_hat)) %>%
  mutate(sigma_lab = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = delta/(sigma^2)) %>%
  mutate(SampSize_lab = paste0("n= ", SampSize))

plot_est <- df %>%
  ggplot() +
  aes(x=Ratio, y = (Variance - sigma^2)/sigma^2, colour = as.factor(SampSize)) +
  geom_line(size = 1.4) +
  geom_hline(aes(yintercept = 0), colour = "#6C0E23", size = 1.4) + 
  scale_color_manual(name = "Sample Size",
                     values = colorRampPalette(c("#caf0f8", "#00b4d8","#03045e"))(length(n)),
                     labels = paste0("n=", n)) +
  facet_wrap(~sigma_lab, labeller=label_parsed) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  xlab(TeX(r'(Ratio $\delta/\sigma^2$)')) +
  ylab(TeX(r'(${(\sigma^2 - \hat{\sigma^2})}/{sigma^2}$)')) +
  theme_classic() +
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24), 
        strip.text = element_text(size = 24),
        plot.tag = element_text(face = "bold", size = 24)) +
  NULL

ggsave(filename = "Supplementary/figures/VarianceEstimationSampleSize.pdf",
       width = 220, 
       height = 150, 
       units = "mm", 
       dpi = 600)
