# --------------------------------- Figure 1 --------------------------------- #

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


#-- Set paramaters

n <- 300
sd1 <- 1
sd2 <- 2
cov_mat <- list(diag(sd1^2, 2),
                cbind(c(sd2^2, 0), c(0, sd1^2)))
delta <- 10
tau <- 0.4
true_classes <- as.factor(rep(1:2, each = n/2))
alpha <- .05

nsimu <- 1000
set.seed(17012024)
seeds <- runif(nsimu, 1, 10^9)

#-- Colors
cluster_col <- c("#294122", "#EB3D00", "#FFBBA6")
results_col <- c("#334EAC", "#BAD6EB", "#6A2A23")
# Panel A: Data representation 
set.seed(2024)
X <- data.frame(X1=c(rnorm(n/2, mean = 0, sd = sd1),
                     rnorm(n/2, mean = delta, sd = sd2)),
                X2 = c(rnorm(n/2, mean = 0, sd = sd1),
                       rnorm(n/2, mean = 0, sd = sd1)))

cl_X <- km_fun(X, K=3)
plt1 <- ggplot(X) + 
  aes(x=X1, y=X2) + 
  geom_density_2d(aes(colour="Original data"), size = 1.2, alpha = .5) +
  scale_colour_manual(name = "", values = "black") +
  ggnewscale::new_scale_colour() +
  geom_point(aes(colour = cl_X), size = 4) +
  scale_colour_manual(name = "Clusters", 
                      values = cluster_col,
                      labels = c(TeX(r'($C_1$)'),
                                 TeX(r'($C_2$)'),
                                 TeX(r'($C_3)'))) +
  xlab(TeX(r'($X_1$)')) +
  ylab(TeX(r'($X_2$)')) +
  theme_classic() +
  # theme(legend.position = "bottom") +
  NULL

plt1



#Simulations 

res <- pblapply(1:nsimu, function(s){
  set.seed(seeds[s])
  X <- data.frame(X1=c(rnorm(n/2, mean = 0, sd = sd1),
                       rnorm(n/2, mean = delta, sd = sd2)),
                  X2 = c(rnorm(n/2, mean = 0, sd = sd1),
                         rnorm(n/2, mean = 0, sd = sd1)))
  
  cl_X <- km_fun(X, K=3)
  
  # Data Fission with variance using the true classes
  fiss_true <- intra_cluster_fission(X, tau = tau, 
                                     cl_ref = true_classes,
                                     sigma_c = cov_mat)
  cl_fiss_true <- km_fun(fiss_true$fX, K=3)
  
  clToTest_true <- order_cluster(fiss_true$gX[,1], cl_fiss_true)
  
  ttest_res_true <- t.test(fiss_true$gX[cl_fiss_true == clToTest_true[1],1], fiss_true$gX[cl_fiss_true== clToTest_true[2],1])$p.value
  ttest_res_true
  # Data Fission with global estimate variance
  fiss_glob <- data_fission(X = as.matrix(X), tau = tau, Sigma = cov(X))
  cl_fiss_glob <- km_fun(fiss_glob$fX, K=3)
  
  clToTest_glob <- order_cluster(fiss_glob$gX[,1], cl_fiss_glob)
  
  ttest_res_glob <- t.test(fiss_glob$gX[cl_fiss_glob ==clToTest_glob[1]], fiss_true$gX[cl_fiss_glob==clToTest_glob[2]])$p.value
  
  
  # Data Fission with intra-cluster covariance but estimating using a first clustering on X
  fiss_est <- intra_cluster_fission(X, tau = tau, cl_ref = cl_X)
  cl_fiss_est <- km_fun(fiss_est$fX, K=3)
  
  clToTest_est <- order_cluster(fiss_est$gX[,1], cl_fiss_est)
  ttest_res_est <- t.test(fiss_glob$gX[cl_fiss_est ==clToTest_est[1]], fiss_true$gX[cl_fiss_est ==clToTest_est[2]])$p.value
  
  # Output results 
  
  df_out <- data.frame(Pvalues = c(ttest_res_true,
                                   ttest_res_glob,
                                   ttest_res_est),
                       Method = rep(c("True Intra Cluster Variance",
                                      "Global Variance", 
                                      "Estimated intra-cluster variance")))
  return(df_out)
}, cl = 5)

res_data <- do.call(rbind.data.frame, res)
plt2 <- do.call(rbind.data.frame, res) %>% 
  ggplot() + 
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + xlab("Theoretical Quantiles") + 
  stat_qq(aes(sample = Pvalues, colour = factor(Method, levels = c("True Intra Cluster Variance",
                                                                   "Estimated intra-cluster variance",
                                                                   "Global Variance"))), distribution = qunif, size = 3) +
  scale_colour_manual(name = "Variance", 
                      values = results_col,
                      labels = c(TeX(r'($\Sigma_g$)'),
                                 TeX(r'($\hat{\Sigma}_{\hat{g}}$)'),
                                 TeX(r'($\hat{\Sigma}$)'))) +
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + theme_classic() 
plt2

plt_final <- plt1 + plt2 +
  plot_layout(widths = c(10, 12)) + 
  plot_annotation(tag_levels = "A") &
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        plot.tag = element_text(face = "bold"))
plt_final
ggsave(plt_final, filename = "figures/figure1.pdf",
       width = 350, 
       height = 120, 
       units = "mm",
       dpi = 600)
