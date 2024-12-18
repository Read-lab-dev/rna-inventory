library(CountClust)
library(data.table)
library(tools)
library(tidyr)

## preprocess cm for nmf analysis
## @param cm log transformed and centered cm
lda_df_preprocessing <- function(cm){
  nonzero_rows <- rowSums(cm)!=0
  nonzero_cols <- colSums(cm)!=0
  cm <- cm[nonzero_rows, nonzero_cols]
  message("The dimension of current cm is ", dim(cm))
  return(cm)
}

findOptimalK <- function(cm, save_path, n_topics=seq(2,10), tolerance=0.1, num_trials=5){
  ## Prepare dir to store outputs 
  sub_dir <- lapply(n_topics, function(x) paste0(x, "topics_tol", tolerance))
  names(sub_dir) <- paste0("topic_", n_topics)
  models <- list()
  
  ## Fit models using different n_topics
  for (i in seq(length(n_topics))){
    current_dir = paste0(save_path, '/', sub_dir[i])
    if (!dir.exists(current_dir)){
      dir.create(current_dir)
    }
    message("Fit modeling...")
    FitGoM(t(as.matrix(cm)), K=n_topics[i], tol=tolerance, num_trials=num_trials,
           path_rda=paste0(save_path, '/', sub_dir[i], '/FitGoM_k', n_topics[i], '_tol', tolerance, '.rda'))
  }
  
  ## Load models into a list 
  message("Loading models...")
  for (m in n_topics) {
    load(paste0(save_path, '/', sub_dir[[paste0("topic_", m)]], '/FitGoM_k', m, '_tol', tolerance, '.rda'))
    name = paste0(m,"_topics")
    models[[name]] <- Topic_clus 
  }
  
  ## Compute diagnostic statistics 
  message("Running compGoM...")     
  out <- compGoM(t(as.matrix(cm)), models)
  message("Saving output...")
  saveRDS(out, file=paste0(save_path, "/compGoM_output.rds"))
  
  ## calculate and plot AIC and plot AIC and BIC
  names(out) <- names(sub_dir)
  message("Plotting AIC and BIC...")
  # Collect measurements from output 
  bic_plot <- sapply(names(out), function(x) out[[x]]$BIC)
  loglik_plot <- sapply(names(out), function(x) out[[x]]$loglik)
  AIC <- 2*(n_topics-1)-2*loglik_plot
  aicbic <- data.frame(AIC=AIC, BIC=bic_plot, k=n_topics)
  ab <- gather(aicbic, type, value, c(AIC, BIC))
  ggplot(ab, aes(x=k, y=value, color=type)) +
    geom_point(size=3) +
    scale_color_manual(name="measure", values = c("AIC"="orangered", "BIC"="dodgerblue4")) +
    xlab("no. of topics")
  ggsave(paste0(save_path, "/aicbic_plot.png"), width=5, height=4)
}







