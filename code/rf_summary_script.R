rf44.container <- list()
rf44.container[[1]] <- readRDS("~/Research/CMPSpaper/code/saved_rds/h44_rf2_features_rmo.rds")
  readRDS("~/Research/CMPSpaper/code/saved_rds/h44_rf_features_def.rds")
rf44.container[[2]] <- 
  readRDS("~/Research/CMPSpaper/code/saved_rds/h44_rf_legacyfeatures_def.rds")

df.rf <- lapply(rf44.container, function(rf_result) {
  with(rf_result$rf.table[[1]], {
    tibble(rf.diff, rf.diff.med, rf.max, rf.maxbar) %>% 
      apply(2, compute_var_ratio_anova, label = type_truth, MS = FALSE)
  })
}) %>% do.call(rbind, .) %>% as.data.frame()

df.rf$rf <- c("features", "legacy_feature")
df.rf$span3 <- c("default", "default")

rf252.container <- list()
rf252.container[[1]] <- readRDS("~/Research/CMPSpaper/code/saved_rds/h252_rf2_features_rmo.rds")
  readRDS(filepath)
  readRDS("~/Research/CMPSpaper/code/saved_rds/h252_rf_features_def.rds")
rf252.container[[2]] <- 
  readRDS("~/Research/CMPSpaper/code/saved_rds/h252_rf_result_25.rds")
rf252.container[[3]] <- 
  readRDS("~/Research/CMPSpaper/code/saved_rds/h252_rf_result_75.rds")

df.rf.252 <- lapply(rf252.container, function(rf_result) {
  with(rf_result$rf.table[[1]], {
    tibble(rf.diff, rf.diff.med, rf.max, rf.maxbar) %>% 
      apply(2, compute_var_ratio_anova, label = type_truth, MS = FALSE)
  })
}) %>% do.call(rbind, .) %>% as.data.frame()

metric_plot_helper(rf252.container[[1]]$rf.table[[1]], "rf.diff.med", scaled = TRUE)
metric_plot_helper(rf44.container[[1]]$rf.table[[1]], "rf.diff", scaled = TRUE)


df.rf.252$rf <- c("features", "legacy_feature", "legacy_feature")
df.rf.252$span3 <- c("default", "25", "75")

df.rf

df.rf.252
