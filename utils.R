library(dplyr)

km_fun <- function(x, K){
  km <- kmeans(x, centers = K, nstart = 1)
  return(as.factor(km$cluster))
}

intra_cluster_fission <- function(X, tau=1, cl_ref, sigma_c = NULL){
  #cl_ref: A reference clustering used to compute intra-cluster variance
  # sigma_c: a list containing intra-cluster cov if knwon
  if (!is.factor(cl_ref)){
    cl_ref <- as.factor(cl_ref)
  }
  if (is.null(sigma_c)){
    sigma_c <- lapply(levels(cl_ref), function(c){
      cov(X[cl_ref == c,])
    })
  }
  intra_cluster_res <- lapply(levels(cl_ref), function(c){
    
    fission_c <- data_fission(as.matrix(X)[cl_ref==c,], Sigma = sigma_c[[c]], tau = tau)
  })
  
  fX <- lapply(intra_cluster_res, function(l){
    return(l$fX)
  })
  
  gX <- lapply(intra_cluster_res, function(l){
    return(l$gX)
  })
  
  return(list(fX =  do.call('rbind',fX),
              gX =  do.call('rbind',gX)))
}

order_cluster <- function(x, cl){
  # x : the variable to test (where the cluster should be ordered) 
  # cl : A three clusters partitions of the data
  
  df <- data.frame(x = x,
                   Cluster = as.factor(cl))
  
  ord_mean <- df %>% group_by(Cluster) %>% summarise(M = mean(x))
  clDiff <- abs(c(ord_mean$M[1]-ord_mean$M[2],
                  ord_mean$M[1]-ord_mean$M[3],
                  ord_mean$M[2] - ord_mean$M[3]))
  
  cl.ord <- which.min(clDiff)
  if(cl.ord == 1){
    cl1 <- 1
    cl2 <- 2 
  }else if(cl.ord == 2){
    cl1 <- 1
    cl2 <- 3
  }else{
    cl1 <- 2
    cl2 <- 3
  }
  return(c(cl1, cl2))
}
