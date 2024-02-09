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
n <- 500
true_class <- rep(1:2, each = n/2)
nsimu <- 1000

#-- Simulations with the true number of classes 
sim_res <- pblapply(1:nsimu, function(ns){
  set.seed(310123*ns)
  x1 <- rnbinom(n=n/2, prob = 0.5, size = 5)
  x2 <- rnbinom(n=n/2, prob = 0.4, size = 40)
  y <- rnbinom(n=n, prob = .5, size = 5)
  X <- cbind(c(x1,x2), 
             y)

  #--- Data thinning with the global parameter
  #-- Overdisperion estimation 
  overdisp_global <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[,p], mu = mean(X[,p]))
  })
  
  res <- countsplit::countsplit(X, overdisps = overdisp_global)
  clust <- km_fun(as.matrix(res[[1]]), K=2)
  ari <- mclust::adjustedRandIndex(clust, true_class)
  
  #--- Data thinning with the true intra-class parameter 
  #-- Overdisperion estimation 
  
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
  
  clust_intra <- km_fun(Xtrain, K=2)
  ari_intra <- mclust::adjustedRandIndex(clust_intra, true_class)
  return(data.frame(ARI = c(ari, ari_intra),
                    Overdispersion = c("Global", "Intra-comp")))
}, cl = 6)


#-- Simulations with the wrong number of classes 
sim_res_wrong <- pblapply(1:nsimu, function(ns){
  set.seed(310123*ns)
  x1 <- rnbinom(n=n/2, prob = 0.5, size = 5)
  x2 <- rnbinom(n=n/2, prob = 0.4, size = 40)
  y <- rnbinom(n=n, prob = .5, size = 5)
  X <- cbind(c(x1,x2), 
             y)
  
  #--- Data thinning with the global parameter
  #-- Overdisperion estimation 
  overdisp_global <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[,p], mu = mean(X[,p]))
  })
  
  res <- countsplit::countsplit(X, overdisps = overdisp_global)
  clust <- km_fun(as.matrix(res[[1]]), K=3)
  clToTest <- order_cluster(as.matrix(res[[1]])[,1], clust)
  res_test <- wilcox.test(as.matrix(res[[2]])[clust == clToTest[1],1], 
                          as.matrix(res[[2]])[clust == clToTest[2],1])$p.value
  
  #--- Data thinning with the true intra-class parameter 
  #-- Overdisperion estimation 
  
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
  return(data.frame(pval = c(res_test, res_test_intra),
                    Overdispersion = c("Global", "Intra-comp")))
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
  theme(legend.position = "bottom") +
  NULL


plt2 <- do.call(rbind.data.frame, sim_res_wrong) %>% 
  ggplot() + 
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + xlab("Theoretical Quantiles") + 
  stat_qq(aes(sample = pval, colour = factor(Overdispersion)),
              distribution = qunif, size = 3) +
  scale_colour_manual(name = "Overdispersion", 
                      values = results_col[c(3,2)],
                      labels=c("Global",
                               "True intra-classe")) +
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + theme_classic() 
plt2

plt_final <- plt1 + plt2 +
  plot_layout(widths = c(12, 10)) + 
  plot_annotation(tag_levels = "A") &
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        plot.tag = element_text(face = "bold"))
plt_final

ggsave(plt_final, filename = "figures/figure4.pdf",
       width = 400, 
       height = 150, 
       units = "mm",
       dpi = 600)
 



