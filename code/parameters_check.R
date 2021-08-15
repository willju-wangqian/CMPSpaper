library(MASS)
source("CMPSpaper_writing/func_collection.R")

######
# generate results from results_generation.R
check.idx <- 2

CMPS_hamby252_results$cmps.table[[check.idx]] %>% length()

tmp.table <- CMPS_hamby252_results$cmps.table[[check.idx]]

#########################
# Aug. 4
# modify func for median
compute_avgscore_denoise(rep(1:6, 6),
                         rep(1:6, each = 6), 
                         cp1$cmps_score, FUNC = median)


aa <- tmp.table$cmps.table.m[[1]]
aa$land1idx <- aa$land1 %>% str_sub(-1, -1) %>% as.numeric()
aa$land2idx <- aa$land2 %>% str_sub(-1, -1) %>% as.numeric()
# compute the new CMPS metric
## minus the average of non match
avg.phase <- compute_avgscore_denoise(
  land1 = aa$land1idx,
  land2 = aa$land2idx,
  score = aa$cmps_score,
  addNA = TRUE
)
avg.phase

rbenchmark::benchmark(
  "cmps_denoise" = {
    avg.phase <- compute_avgscore_denoise(
      land1 = aa$land1idx,
      land2 = aa$land2idx,
      score = aa$cmps_score,
      addNA = TRUE
    )
  }, replications = 100)

## with the bullet_land_predict function, we cannot
## guarantee to get samesource
aa$same_source <- 
  with(aa,
       bullet_to_land_predict(
         land1 = land1idx,
         land2 = land2idx,
         cmps_score,
         difference = 1
       ))

with(aa, {
  mean(cmps_score[same_source], na.rm=TRUE) - 
    mean(cmps_score[!same_source], na.rm = TRUE)
})

# write a function to compute the variance ratio
## check the definition in Joe's paper
## 
tab1 <- CMPS_hamby252_results$cmps.table[[check.idx]] %>% 
  dplyr::select(cmps.maxbar, type_truth)
tab1
tab1$cmps.maxbar.scaled <- with(tab1, {
  x <- cmps.maxbar
  (x - mean(x)) / sd(x)
})


# variance ratio
with(tab1, compute_var_ratio(cmps.maxbar, type_truth))
with(tab1, compute_var_ratio(cmps.maxbar.scaled, type_truth))

# linear discriminant analysis
model <- lda(type_truth ~ cmps.maxbar, data = tab1)
model$svd

# anova
mm <- lm(cmps.maxbar ~ type_truth, data = tab1)
ss.aov <- anova(mm)
ss.aov
ss.aov$`Sum Sq`[1] / ss.aov$`Sum Sq`[2]

###########
# compute new metric and compare, using variance ratio

for (i in 1:length(CMPS_hamby252_results$signame)) {
  CMPS_hamby252_results$cmps.table[[i]] -> hamby252.cmps
  
  hamby252.cmps <- hamby252.cmps %>% mutate(
    cmps.denoise = cmps.table.m %>% purrr::map_dbl(
      .f = function(t) {
        landidx1 <-
          t$land1 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        landidx2 <-
          t$land2 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        compute_avgscore_denoise(landidx1, landidx2, t$cmps_score, addNA = TRUE)
      }
    )
  )
  
  CMPS_hamby252_results$cmps.table[[i]] <- hamby252.cmps
}


##############
# regenerate plots
for (i in 1:4) {
  hamby252.cmps <- CMPS_hamby252_results$cmps.table[[i]]
  
  hamby252.plot.list <- list()
  titlee <- CMPS_hamby252_results$titlee[[i]]
  
  hamby252.plot.list[[1]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.max.m,
                       fill = as.factor(type_truth)), binwidth = 1) +
    labs(
      fill = "Comparison Type",
      x = expression(CMPS[max]),
      # title = expression(paste("Hamby252 - ", CMPS[max], " Distribution")),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = seq(0, 27, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    font("x.text", size = 6) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  hamby252.plot.list[[2]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.maxbar.m,
                       fill = as.factor(type_truth)), binwidth = 1) +
    labs(
      x = expression(bar(CMPS)[max]),
      fill = "Comparison Type",
      # title = expression(paste(
      #   "Hamby252 - ", bar(CMPS)[max], " Distribution"
      # )),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = seq(0, 24, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  hamby252.plot.list[[3]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.denoise,
                       fill = as.factor(type_truth)), binwidth = 1) +
    labs(
      x = expression(bar(CMPS)[denoise]),
      fill = "Comparison Type",
      # title = expression(paste(
      #   "Hamby252 - ", bar(CMPS)[max], " Distribution"
      # )),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = seq(0, 24, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  plot <- ggarrange(plotlist = hamby252.plot.list,
                    nrow = 1,
                    ncol = 3,
                    common.legend = TRUE, legend = "bottom")
  plot <- annotate_figure(plot, 
                          top = text_grob(com.title252))
  CMPS_hamby252_results$plot[[i]] <- plot
}
CMPS_hamby252_results$plot[[1]]

# compute variance ratio for the three plots
with(CMPS_hamby252_results$cmps.table[[1]], {
  cbind(cmps.max.m, cmps.maxbar.m, cmps.denoise) %>% 
    apply(2, compute_var_ratio, label = type_truth)
}) 

CMPS_hamby252_results$cmps.table[[1]] %>% filter(type_truth == "KM") %>% 
  dplyr::select(cmps.maxbar.m, cmps.denoise)

#######
# cmps_score / nseg

for (i in 1:length(CMPS_hamby252_results$signame)) {
  CMPS_hamby252_results$cmps.table[[i]] -> hamby252.cmps
  
  hamby252.cmps <- hamby252.cmps %>% mutate(
    cmps.table.m = cmps.table.m %>% purrr::map(
      .f = function(t) {
        t$cmps_score_scaled <- t$cmps_score / t$cmps_nseg
        
        t
      }
    )
  )
  
  hamby252.cmps <- hamby252.cmps %>% mutate(
    cmps.denoise_scaled = cmps.table.m %>% purrr::map_dbl(
      .f = function(t) {
        landidx1 <-
          t$land1 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        landidx2 <-
          t$land2 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        compute_avgscore_denoise(landidx1, landidx2, t$cmps_score_scaled, addNA = TRUE)
      }
    ),
    cmps.max.m_scaled = cmps.table.m %>% purrr::map_dbl(
      .f = function(t) {
        max(t$cmps_score_scaled, na.rm = TRUE)
      }
    ),
    cmps.maxbar.m_scaled = cmps.table.m %>% purrr::map_dbl(
      .f = function(t) {
        landidx1 <-
          t$land1 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        landidx2 <-
          t$land2 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        compute_average_scores(landidx1, landidx2, t$cmps_score_scaled, addNA = TRUE) %>% max()
      }
    )
  )
  
  CMPS_hamby252_results$cmps.table[[i]] <- hamby252.cmps
}

CMPS_hamby252_results$plot_scaled <- list()

for (i in 1:4) {
  hamby252.cmps <- CMPS_hamby252_results$cmps.table[[i]]
  
  hamby252.plot.list <- list()
  titlee <- CMPS_hamby252_results$titlee[[i]]
  
  hamby252.plot.list[[1]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.max.m_scaled,
                       fill = as.factor(type_truth)), binwidth = 0.04) +
    labs(
      fill = "Comparison Type",
      x = expression(CMPS[max]),
      # title = expression(paste("Hamby252 - ", CMPS[max], " Distribution")),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.05)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    font("x.text", size = 6) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  hamby252.plot.list[[2]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.maxbar.m_scaled,
                       fill = as.factor(type_truth)), binwidth = 0.04) +
    labs(
      x = expression(bar(CMPS)[max]),
      fill = "Comparison Type",
      # title = expression(paste(
      #   "Hamby252 - ", bar(CMPS)[max], " Distribution"
      # )),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.05)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  hamby252.plot.list[[3]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.denoise_scaled,
                       fill = as.factor(type_truth)), binwidth = 0.04) +
    labs(
      x = expression(bar(CMPS)[denoise]),
      fill = "Comparison Type",
      # title = expression(paste(
      #   "Hamby252 - ", bar(CMPS)[max], " Distribution"
      # )),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.05)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  plot <- ggarrange(plotlist = hamby252.plot.list,
                    nrow = 1,
                    ncol = 3,
                    common.legend = TRUE, legend = "bottom")
  plot <- annotate_figure(plot, 
                          top = text_grob(com.title252))
  CMPS_hamby252_results$plot_scaled[[i]] <- plot
}

data(bullets)
land2_3 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
land1_2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]
cmps_with_multi_scale <- 
  extract_feature_cmps(land2_3$sig, land1_2$sig, include = "full_result",
                       seg_length = 50, npeaks.set = c(10,6,4,2,1))
cmps_with_multi_scale$CMPS.score
cmps_with_multi_scale$nseg


check_idx <- 1

CMPS_hamby252_results$plot_scaled[[check_idx]]
CMPS_hamby252_results$plot[[check_idx]]
CMPS_hamby252_results$titlee[[check_idx]]
# compute variance ratio for the three plots
with(CMPS_hamby252_results$cmps.table[[check_idx]], {
  cbind(cmps.max.m, cmps.maxbar.m, cmps.denoise, 
        cmps.max.m_scaled, cmps.maxbar.m_scaled, cmps.denoise_scaled,
        cmps.max, cmps.maxbar) %>% 
    apply(2, compute_var_ratio, label = type_truth)
}) 



# ANOVA?
# https://en.wikipedia.org/wiki/F-test
with(CMPS_hamby252_results$cmps.table[[check_idx]], {
  cbind(cmps.max.m, cmps.maxbar.m, cmps.denoise, 
        cmps.max.m_scaled, cmps.maxbar.m_scaled, cmps.denoise_scaled) %>% 
    apply(2, function(ss) {
      mm <- lm(ss ~ type_truth)
      ss.aov <- anova(mm)
      ss.aov$`Sum Sq`[1] / ss.aov$`Sum Sq`[2]
    })
}) 

with(CMPS_hamby252_results$cmps.table[[check_idx]], {
  cbind(cmps.max.m, cmps.maxbar.m, cmps.denoise, 
        cmps.max.m_scaled, cmps.maxbar.m_scaled, cmps.denoise_scaled) %>% 
    apply(2, function(ss) {
      mm <- lm(ss ~ type_truth)
      ss.aov <- anova(mm)
      ss.aov$`Mean Sq`[1] / ss.aov$`Mean Sq`[2]
    })
}) 

#
with(CMPS_hamby252_results$cmps.table[[check_idx]], {
  cbind(cmps.max.m, cmps.maxbar.m, cmps.denoise, 
        cmps.max.m_scaled, cmps.maxbar.m_scaled, cmps.denoise_scaled) %>% 
    apply(2, compute_var_ratio_anova, label = type_truth, MS = TRUE)
}) 

with(CMPS_hamby252_results$cmps.table[[check_idx]], {
  cbind(cmps.max.m, cmps.maxbar.m, cmps.denoise, 
        cmps.max.m_scaled, cmps.maxbar.m_scaled, cmps.denoise_scaled) %>% 
    apply(2, compute_var_ratio_anova, label = type_truth, MS = FALSE)
}) 


















