library(Seurat)
library(DropletUtils)
library(patchwork)
library(zellkonverter)
library(scuttle)
library(scater)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(FactoMineR)
library(factoextra)
library(umap)

theme_set(theme_bw())

#-- Data preprocessing 
#--- Load data
set.seed(2024)
TS_bone_marrow <- zellkonverter::readH5AD("~/Downloads/TS_Bone_Marrow.h5ad",
                                          verbose = T,
                                          layers = F, varm = F, obsm = F, varp = F, obsp = F, uns = F
)
names(assays(TS_bone_marrow)) <- "counts"

#--- Log-normalisation 
log_norm <- logcounts(logNormCounts(TS_bone_marrow))

#--- Selection of the cells of interests 
counts_select_seurat <- TS_bone_marrow[, which((TS_bone_marrow$donor == "TSP14") & 
                                          (TS_bone_marrow$cell_ontology_class %in% c("memory b cell", 
                                                                                     "monocyte",
                                                                                     "nk cell",
                                                                                     "macrophage")))]


counts_select <- t(counts(counts_select_seurat))
dim(counts_select)

#--- Log-normalisation
log_norm_select <- log_norm[,which((TS_bone_marrow$donor == "TSP14") & 
                                      (TS_bone_marrow$cell_ontology_class %in% c("memory b cell", 
                                                                                 "monocyte",
                                                                                 "nk cell",
                                                                                 "macrophage")))]
log_counts <- t(log_norm_select)
stopifnot(dim(log_counts) == dim(counts_select))

#--- Feature selection based on highly variable genes 
var_counts <- apply(log_counts, 2, var)
mean_counts <- apply(log_counts, 2, mean)
X <- log_counts[, order(var_counts, decreasing = T)[1:500]]

X_raw <- counts_select[, which(colnames(counts_select) %in% colnames(X))]
stopifnot(dim(X_raw) == dim(X))

#--- Data description
cell_lab <- TS_bone_marrow$cell_ontology_class[which(
  (TS_bone_marrow$donor == "TSP14") &
    (TS_bone_marrow$cell_ontology_class %in%
       c("memory b cell", "monocyte", "nk cell", "macrophage")))]
stopifnot(length(cell_lab) == nrow(X))
cell_col <- c('#5C0029', 
              "#E63946",
              "#A8DADC",
              "#457B9D")
levels(cell_col) <- unique(cell_lab)
pca_all <- PCA(X, graph = F, ncp = 50, scale.unit = T)
fviz_pca_ind(pca_all, geom = "point", col.ind = cell_lab, axes = c(1,2))

pca_all$ind$coord %>% data.frame() %>%
  mutate(CellPop = cell_lab) %>%
  ggplot() +
  aes(x=Dim.1, y = Dim.2, colour = CellPop) +
  geom_point(size = 2, alpha = .8) +
  geom_hline(yintercept = 0, size = .9, colour = "black", linetype = 2) +
  geom_vline(xintercept = 0, size = .9, colour = "black", linetype = 2) +
  scale_colour_manual(name = "Cell Population",
                      values =  cell_col) +
  xlab(paste0("Dim 1 (", round(pca_all$eig[1,2],2), "%)")) +
  ylab(paste0("Dim 2 (", round(pca_all$eig[2,2],2), "%)")) +
  theme_classic() +
  NULL

#-- Overdispersion estimate 

#--- Global

overdisp_global_npreg <- sapply(1:ncol(X_raw), function(p){
  res <- try(npreg::theta.mle(X_raw[,p], mu = mean(X_raw[,p])))
  if(class(res) == "try-error"){
    res_temp <- NA
  }
  else{
    res_temp <- res
  }
  return(res_temp)
})

overdisp_global_vst <- sctransform::vst(t(X_raw)) 

plot_aggrement <- data.frame(npreg = overdisp_global_npreg,
                            VST = overdisp_global_vst$model_pars[,1]) %>%
  ggplot() +
  aes(x=npreg, y = VST) +
  geom_point(size = 2, alpha = .5) + 
  xlab("Overdispersion estimate with npreg") +
  ylab("Overdispersion estimate with VST") +
  # scale_x_log10() +
  # scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, colour = "darkred", size = 1.2) +
  NULL  

plot_aggrement_zoom <- data.frame(npreg = overdisp_global_npreg,
                                  VST = overdisp_global_vst$model_pars[,1]) %>%
  ggplot() +
  aes(x=npreg, y = VST) +
  geom_point(size = 2, alpha = .5) + 
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  xlab("Overdispersion estimate with npreg") +
  ylab("Overdispersion estimate with VST") +
  geom_abline(slope = 1, intercept = 0, colour = "darkred", size = 1.2) +
  NULL  

plot_aggrement + plot_aggrement_zoom 


#--- Intra monocyte
compute_overdisp <- function(cell_pop){
  overdisp_intra <- pbapply::pbsapply(1:ncol(X_raw), function(p){
    res <- try(npreg::theta.mle(X_raw[cell_lab == cell_pop,p], mu = mean(X_raw[cell_lab == cell_pop,p])))
    
    if(class(res) == "try-error"){
      res_temp <- NA
    }
    else{
      res_temp <- res
    }
    return(res_temp)
  }, cl = 5)
  df_overdisp_intra <- data.frame(Gene = colnames(X_raw),
                                       Intra = overdisp_intra)
  
  df_overdisp_glob <- data.frame(Gene = colnames(X_raw),
                                 Global = overdisp_global_npreg)
  
  df_all <- merge(df_overdisp_glob, df_overdisp_intra, by = "Gene")
  df_all$CellType <-  cell_pop
  
  plt <- merge(df_overdisp_intra, df_overdisp_glob, by = "Gene") %>%
    ggplot() + 
    aes(x=Global, y = Intra) +
    geom_point(size = 2, alpha = .5, colour = cell_col[which(levels(cell_col)==cell_pop)]) + 
    xlim(c(0, 4)) +
    ylim(c(0, 4)) +
    xlab(paste0("Intra-", cell_pop, " overdispersion ")) +
    ylab("Global Overdispersion") +
    geom_abline(slope = 1, intercept = 0, colour = "darkred", size = 1.2) +
    NULL  
  R2 <- merge(df_overdisp_glob, df_overdisp_intra, by = "Gene") %>% filter(Global < 4 & Intra < 4) %>% 
    summarise(R2 = cor(Global, Intra))
  return(list(plot = plt, R2 = R2, overdisp = df_all))
}

res <- lapply(unique(cell_lab), function(c){compute_overdisp(c)})

overdispersion_estimate <- do.call("rbind.data.frame",
                                   lapply(1:length(res), function(l){res[[l]]$overdisp}))
  
write.csv(overdispersion_estimate, 
          file = "results/overdispersion_estimate.csv",
          row.names = F)
all_plt <- lapply(1:length(res), function(l){res[[l]]$plot})
do.call(cowplot::plot_grid, all_plt)

all_R2 <- lapply(1:length(res), function(l){res[[l]]$R2})
do.call(rbind.data.frame, all_R2) %>% mutate(`Cell population` = unique(cell_lab))
