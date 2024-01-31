# --------------------------------- Figure 3 --------------------------------- #

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


sim_fun <- function(seed, delta=20, n=100, sd = 1, tau = 0.4){
  set.seed(250124*seed)
  X <- c(rnorm(n/2, mean = 0, sd = sd),
         rnorm(n/2, mean = delta, sd = sd))
  variance <- local_var(X)
  Z <- sapply(variance, function(s){rnorm(1, mean = 0, sd = sqrt(s))})
  fX  <- X + tau*Z 
  gX <- X - (1/tau)*Z 
  
  cl <- kmeans(fX, centers = 3, nstart = 100)$cluster
  df <- data.frame(gX = gX, 
                   Cluster = as.factor(cl))
  
  clToTest <- order_cluster(gX, cl)
  res.ttest <- t.test(gX[cl ==clToTest[1]], gX[cl==clToTest[2]])$p.value
  return(data.frame(pval = res.ttest,
                    Variance_hat = mean(variance)))
}

apply_sim_delta <- function(seed, delta_grid, n=100, sd = 1, tau = 0.4){
  res_delta <- lapply(delta_grid, function(d){
    res_temp <- sim_fun(seed = seed, delta = d, sd = sd)
    return(res_temp)})
  temp <- do.call("rbind.data.frame", res_delta)
  temp$delta <- delta_grid
  return(temp)
}

apply_sim_delta_sigma <- function(seed, delta_grid, n=100, sd_grid, tau = 0.4){
  res_sigma <- lapply(sd_grid, function(s){
    res_temp <- apply_sim_delta(seed = seed, delta_grid = delta_grid, n = n, sd = s, tau = .4)
    res_temp$sigma <- s
    return(res_temp)
  })
  temp <- do.call("rbind.data.frame", res_sigma)
  return(temp)
}

#--Param 
n <- 100
delta_grid <- c(seq(0, 3, length.out = 25), seq(3.5, 100, length.out = 25))
tau <- .4
sigma <- c(0.1, 0.5, 1, 2) 
nsimu <- 1000

res_sim <- pblapply(1:nsimu, function(ns){
  res <- apply_sim_delta_sigma(seed = ns, delta_grid = delta_grid, n = n, sd_grid = sigma, tau = tau)
  return(res)
  }, cl = 5)

df <- do.call(rbind.data.frame,res_sim) %>%
  group_by(delta, sigma) %>%
  summarise(TypeI = mean(pval < 0.05), 
            Variance= mean(Variance_hat),
            sdVariance_hat = sd(Variance_hat)) %>%
  mutate(sigma_lab = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = delta/(sigma^2))

plot_typeI <- df %>% 
  ggplot() + 
  aes(x=Ratio, y = TypeI, colour = sigma_lab) +
  # geom_point(size = 3) +
  geom_line(size = 1.4) +
  scale_colour_manual(name = TeX(r'($\sigma^2$)'),
                      values = c("#93B5C6", "#DBC2CF","#998bc0", "#BD4F6C"),
                      labels = c(0.01, 0.5, 1, 4)) +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal levels"),
             linetype = 2,
             size = 1.2) +
  scale_colour_manual(name = "",
                      values = "#6C0E23") +
  scale_x_log10() +
  xlab(TeX(r'(Ratio $\delta/\sigma^2$)')) +
  annotation_logticks(sides = "b") +
  ylab("Empirical Type I error rate") +
  theme_classic() +
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24), 
        strip.text = element_text(size = 24),
        plot.tag = element_text(face = "bold", size = 24)) +
  NULL

plot_est <- df %>% mutate(Upper = sigma^2 + (Variance + sdVariance_hat), 
                          Lower = sigma^2 - (Variance - sdVariance_hat)) %>%
  ggplot() +
  aes(x=Ratio, y = (Variance - sigma^2)/sigma^2, colour = sigma_lab) +
  geom_line(size = 1.4) +
  scale_colour_manual(name = TeX(r'($\sigma^2$)'),
                     values = c("#93B5C6", "#DBC2CF", "#998bc0", "#BD4F6C"),
                     labels = c(0.01, 0.5, 1, 4)) +
  geom_hline(aes(yintercept = 0), colour = "#6C0E23", size = 1.4) + 
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

plot_est/plot_typeI + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 24))
ggsave(filename = "figures/figure3.pdf", 
       width = 220,
       height = 175, 
       units = "mm",
       dpi = 600)
 