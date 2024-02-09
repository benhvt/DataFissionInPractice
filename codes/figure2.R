# --------------------------------- Figure 2 --------------------------------- #

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
sim_fun <- function(n,sigma, sigma_hat, tau, nsim){
  res <- pbsapply(1:nsim, function(ns){
    X <- rnorm(n, mean = 0, sd = sigma)
    Z <- rnorm(n, mean = 0, sd = sigma_hat)
    X1 <- X + tau*Z
    X2 <- X - (1/tau)*Z
    
    cl <- km_fun(X1, K=2)
    test <- t.test(X2~cl)
    return(test$p.value)
  }, cl = parallel::detectCores()-1)
  return(mean(res < 0.05))
}

sim_fun_sigma_grid <- function(n,sigma, sigma_hat, tau, nsim){
  res <- pbsapply(sigma_hat, function(s){sim_fun(n, sigma, sigma_hat = s, tau, nsim)})
  return(res)
}

sim_fun_sigma <- function(n,sigma, sigma_hat, tau, nsim){
  res <- pblapply(sigma, function(s){sim_fun_sigma_grid(n, sigma = s, sigma_hat, tau, nsim)})
  return(res)
}

#-- Simulations 
n <- 100
sigma <- c(0.1, 0.5, 1, 2)
sigma_grid <- sort(c(seq(0, 4, length.out = 50), 2, 0.1))
nsimu <- 1000
tau <- .4

sim_res <- sim_fun_sigma(n=n, sigma = sigma, sigma_hat = sigma_grid, tau = tau, nsim = nsimu) 

#-- Plot 
data.frame(EmpTypeI = unlist(sim_res),
           sigma = rep(sigma, each = length(sigma_grid)),
           sigma_hat = rep(sigma_grid, length(sigma))) %>%
  mutate(TheTypeI = compute_typeI(n = n, tau = tau, sigma = sigma, sigma_hat = sigma_hat, alpha = .05)) %>%
  mutate(sigma_name = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = (sigma^2-sigma_hat^2)/sigma^2) %>%
  ggplot() + aes(x=Ratio, y = EmpTypeI) +
  geom_point(shape = 2, size = 3) +  
  geom_line(aes(x=Ratio, y = TheTypeI, colour = "Theoritical Type I error rate"), linewidth = 1.5) +
  scale_colour_manual(name = " ",
                      values ="#5C7AFF") +
  xlab(TeX(r'($(\sigma^2 - \hat{\sigma^2})/\sigma^2$)')) +
  ylab("Type I error rate") +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal levels"),
             linetype = 2,
             size = 1.2) +
  scale_colour_manual(name = "",
                      values = "#6C0E23") +
  facet_wrap(~as.factor(sigma_name), 
             labeller = label_parsed,
             scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24), 
        strip.text = element_text(size = 24),
        plot.tag = element_text(face = "bold", size = 24)) +
  NULL

ggsave(filename = "figures/figure2.pdf",
       width = 200, 
       height = 160, 
       units = "mm",
       dpi = 600)

##- Second prorposition 
data.frame(EmpTypeI = unlist(sim_res),
           sigma = rep(sigma, each = length(sigma_grid)),
           sigma_hat = rep(sigma_grid, length(sigma))) %>%
  mutate(TheTypeI = compute_typeI(n = n, tau = tau, sigma = sigma, sigma_hat = sigma_hat, alpha = .05)) %>%
  mutate(sigma_name = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = (sigma^2-sigma_hat^2)/sigma^2) %>%
  ggplot() + aes(x=Ratio, y = EmpTypeI, colour = sigma_name) +
  geom_point(shape = 2, size = 4) +  
  scale_colour_manual(name = 'Empirical Type I error rate',
                      values = c("#93B5C6", "#DBC2CF", "#998bc0", "#BD4F6C"),
                      labels = c(TeX(r'($\sigma^2 = 0.01$)'),
                                          TeX(r'($\sigma^2 = 0.25$)'),
                                          TeX(r'($\sigma^2 = 1$)'),
                                          TeX(r'($\sigma^2 = 4$)'))) +
  ggnewscale::new_scale_colour() +
  geom_line(aes(x=Ratio, y = TheTypeI, colour = sigma_name), linewidth = 1.5, alpha = .8) +
  scale_colour_manual(name = "Theoritical Type I error rate",
                      values = c("#677f8b", "#998891", "#6b6186", "#84374c"),
                      labels = c(TeX(r'($\sigma^2 = 0.01$)'),
                                 TeX(r'($\sigma^2 = 0.25$)'),
                                 TeX(r'($\sigma^2 = 1$)'),
                                 TeX(r'($\sigma^2 = 4$)'))) +
  xlab(TeX(r'($(\sigma^2 - \hat{\sigma^2})/\sigma^2$)')) +
  xlim(c(-5,2)) +
  ylab("Type I error rate") +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal levels"),
             linetype = 2,
             size = 1.2) +
  scale_colour_manual(name = "",
                      values = "#6C0E23") +
  NULL

ggsave(filename = "figures/figure2Option2.pdf",
       width = 180, 
       height = 100, 
       units = "mm",
       dpi = 600)
 