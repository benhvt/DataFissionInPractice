# --------------------------------- Figure 4 --------------------------------- #
# --------------------- Negative Binomial Simulations -------------------------# 
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

#-- Parameters 
n <- 100
true_class <- rep(1:2, each = n/2)
nsimu <- 1000

#-- Simulations with the wrong number of classes 
sim_res_wrong <- pblapply(1:nsimu, function(ns){
  set.seed(310123*ns)
  x1 <- rnbinom(n=n/2, prob = 0.5, size = 5)
  x2 <- rnbinom(n=n/2, prob = 0.4, size = 40)
  y <- rnbinom(n=n, prob = .5, size = 5)
  X <- cbind(c(x1,x2), 
             y)
  
  #--- Data thinning with the global parameter
  overdisp_global <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[,p], mu = mean(X[,p]))
  })
  
  res <- countsplit::countsplit(X, overdisps = overdisp_global)
  clust <- km_fun(as.matrix(res[[1]]), K=3)
  clToTest <- order_cluster(as.matrix(res[[1]])[,1], clust)
  res_test <- wilcox.test(as.matrix(res[[2]])[clust == clToTest[1],1], 
                          as.matrix(res[[2]])[clust == clToTest[2],1])$p.value
  
  #--- Data thinning with the true intra-class parameter 
  overdisp_C1 <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[true_class==1,p], mu = mean(X[true_class==1,p]))
  })
  
  overdisp_C2 <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[true_class==2,p], mu = mean(X[true_class==2,p]))
  })
  overdisp_intra <- list(overdisp_C1, overdisp_C2)
  res_temp <- lapply(1:2, function(cl){
    res_datathin <- countsplit::countsplit(X[true_class==cl, ], overdisps = overdisp_intra[[cl]])
  })
  
  Xtrain <- rbind(as.matrix(res_temp[[1]][[1]]),
                  as.matrix(res_temp[[2]][[1]]))
  Xtest <- rbind(as.matrix(res_temp[[1]][[2]]),
                 as.matrix(res_temp[[2]][[2]]))
  
  clust_intra <- km_fun(Xtrain, K=3)
  clToTest_intra <- order_cluster(Xtrain[,1], clust_intra)
  res_test_intra <- wilcox.test(Xtest[clust_intra==clToTest_intra[1],1],
                          Xtest[clust_intra==clToTest_intra[2],1])$p.value  
  
  #--- Data thinning with the wrong intra-cluster parameter 
  cl_ref <- km_fun(X, K=3)
  
  overdisp_C1_hat <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[cl_ref==1,p], mu = mean(X[cl_ref==1,p]))
  })
  
  overdisp_C2_hat <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[cl_ref==2,p], mu = mean(X[cl_ref==2,p]))
  })
  
  overdisp_C3_hat <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[cl_ref==3,p], mu = mean(X[cl_ref==3,p]))
  })
  
  overdisp_intra_hat <- list(overdisp_C1_hat, overdisp_C2_hat, overdisp_C3_hat)
  res_temp_hat <- lapply(1:3, function(cl){
    res_datathin <- countsplit::countsplit(X[cl_ref==cl, ], overdisps = overdisp_intra_hat[[cl]])
  })
  
  Xtrain_hat <- rbind(as.matrix(res_temp_hat[[1]][[1]]),
                  as.matrix(res_temp_hat[[2]][[1]]),
                  as.matrix(res_temp_hat[[3]][[1]]))
  Xtest_hat <- rbind(as.matrix(res_temp_hat[[1]][[2]]),
                 as.matrix(res_temp_hat[[2]][[2]]),
                 as.matrix(res_temp_hat[[3]][[2]]))
  
  
  clust_intra_hat <- km_fun(Xtrain_hat, K=3)
  clToTest_intra_hat <- order_cluster(Xtrain_hat[,1], clust_intra_hat)
  res_test_intra_hat <- wilcox.test(Xtest_hat[clust_intra_hat==clToTest_intra_hat[1],1],
                                Xtest_hat[clust_intra_hat==clToTest_intra_hat[2],1])$p.value  
  
  return(data.frame(pval = c(res_test, res_test_intra, res_test_intra_hat),
                    Overdispersion = c("Global", "Intra-comp", "Intra-cluster")))
}, cl = 6)


#-- Make figure
cluster_col <- c("#294122", "#EB3D00", "#FFBBA6")
results_col <- c("#334EAC", "#BAD6EB", "#6A2A23")
set.seed(09022024)

x1 <- rnbinom(n=n/2, prob = 0.5, size = 5)
x2 <- rnbinom(n=n/2, prob = 0.4, size = 40)
y <- rnbinom(n=n, prob = .5, size = 5)
X <- cbind.data.frame(X1=c(x1,x2), 
           X2=y, 
           TrueClasses=as.factor(true_class))
cl_X <- km_fun(cbind(c(x1,x2), y), K = 3)

plt1 <- ggplot(X) + aes(x=X1, y=X2) + 
  geom_density_2d(aes(colour=TrueClasses), size = 1.2, alpha = .5) +
  scale_colour_manual(name = "True classes", values = c("#274060", "#E70E02")) +
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


plt2 <- do.call(rbind.data.frame, sim_res_wrong) %>% 
  ggplot() + 
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + xlab("Theoretical Quantiles") + 
  stat_qq(aes(sample = pval, colour = factor(Overdispersion,
                                             levels = c('Intra-comp', 'Intra-cluster', 'Global'))),
              distribution = qunif, size = 3) +
  scale_colour_manual(name = "Overdispersion", 
                      values = results_col,
                      labels=c(TeX(r'($\hat{\theta}_{g}$)'),
                               TeX(r'($\hat{\theta}_{\hat{g}}$)'),
                               TeX(r'($\hat{\theta}$)'))) +
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + theme_classic() 
plt2

plt_final <- plt1 + plt2 +
  plot_layout(widths = c(10, 8)) +
  plot_annotation(tag_levels = "A") &
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        # legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        plot.tag = element_text(face = "bold"))
plt_final

ggsave(plt_final, filename = "figures/figure4.pdf",
       width = 350, 
       height = 100, 
       units = "mm",
        dpi = 600)
       



