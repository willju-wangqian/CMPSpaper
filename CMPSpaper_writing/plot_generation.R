library(CMPS)
library(tidyverse)
library(bulletxtrctr)
library(x3ptools)

data(bullets)

lands <- unique(bullets$bulletland)

comparisons <- data.frame(expand.grid(land1 = lands[1:6], land2 = lands[7:12]), 
                          stringsAsFactors = FALSE)

comparisons <- comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, 
                        .f = function(xx, yy) {
                          land1 <- bullets$sigs[bullets$bulletland == xx][[1]]
                          land2 <- bullets$sigs[bullets$bulletland == yy][[1]]
                          land1$bullet <- "first-land"
                          land2$bullet <- "second-land"
                          
                          sig_align(land1$sig, land2$sig)
                        }))

comparisons <- comparisons %>% 
  mutate(cmps = aligned %>% purrr::map(.f = function(a) {
    extract_feature_cmps(a$lands$sig1, a$lands$sig2, include = "full")
  }))


comparisons <- comparisons %>% 
  mutate(
    cmps_score = sapply(comparisons$cmps, function(x) x$CMPS.score),
    cmps_nseg = sapply(comparisons$cmps, function(x) x$nseg)
  )

# dframe <- comparisons %>% select(land1, land2, cmps_score)
# 
# dframe %>% 
#   ggplot(aes(x = land1, y = land2, fill = cmps_score)) +
#   geom_tile() +
#   scale_fill_gradientn(
#     colors=c("darkgrey","darkorange")
#   ) +
#   geom_text(aes(label=cmps_score)) +
#   theme_bw() + 
#   coord_flip() 

# low="black", high = "darkorange"

comparisons <- comparisons %>%
  mutate(
    bulletA = gsub("(\\d)-\\d", "\\1", land1),
    landA = gsub("\\d-(\\d)", "\\1", land1),
    bulletB = gsub("(\\d)-\\d", "\\1", land2),
    landB = gsub("\\d-(\\d)", "\\1", land2)
  )

comparisons %>%   
  group_by(bulletA, bulletB) %>% tidyr::nest() %>%
  mutate(
    cmps_max_bar = data %>% purrr::map_dbl(
      .f = function(d) max(compute_average_scores(land1 = d$landA,
                                                  land2 = d$landB,
                                                  d$cmps_score)))
  )

dframe <- comparisons %>% select(-aligned, -cmps)

dframe %>% ggplot(aes(x = landB, y = landA, fill = cmps_score)) + geom_tile() + 
  scale_fill_gradient2(low = "gray80", high = "darkorange", midpoint = 6) + 
  facet_grid(bulletA ~ bulletB, labeller = "label_both") + xlab("Land B") + 
  ylab("Land A") + theme(aspect.ratio = 1) +
  geom_text(aes(label=cmps_score)) +
  theme_bw()



#######################
# check 1-4 vs 2-5, segment 8 -> expected NA, got INF

#########################
# generate plots for show & tell
CMPS_hamby252_results <- readRDS("~/Research/CMPSpaper/CMPSpaper_writing/data/h252-seg_length.rds")
com.title252 <- expression(paste(
  "Hamby 252 - ", CMPS[max], ", ", bar(CMPS)[max], " and ", bar(CMPS)[diff] ," Distribution"
))


for (i in 1:length(CMPS_hamby252_results$signame)) {
  CMPS_hamby252_results$cmps.table[[i]] -> hamby252.cmps
  
  hamby252.cmps <- hamby252.cmps %>% mutate(
    cmps.max_scaled = cmps.table.m %>% purrr::map_dbl(
      .f = function(t) {
        max(t$cmps_score_scaled, na.rm = TRUE)
      }
    ),
    cmps.maxbar_scaled = cmps.table.m %>% purrr::map_dbl(
      .f = function(t) {
        landidx1 <-
          t$land1 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        landidx2 <-
          t$land2 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        compute_average_scores(landidx1, landidx2, t$cmps_score_scaled, addNA = TRUE) %>% max()
      }
    ),
    cmps.max = cmps.table.m %>% purrr::map_dbl(
      .f = function(t) {
        max(t$cmps_score, na.rm = TRUE)
      }
    ),
    cmps.maxbar = cmps.table.m %>% purrr::map_dbl(
      .f = function(t) {
        landidx1 <-
          t$land1 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        landidx2 <-
          t$land2 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        compute_average_scores(landidx1, landidx2, t$cmps_score, addNA = TRUE) %>% max()
      }
    )
  )
  
  CMPS_hamby252_results$cmps.table[[i]] <- hamby252.cmps
}

# plot
CMPS_hamby252_results$plot_scaled <- list()
CMPS_hamby252_results$plot <- list()

for (i in 1:length(CMPS_hamby252_results$signame)) {
  hamby252.cmps <- CMPS_hamby252_results$cmps.table[[i]]
  
  hamby252.plot.list <- list()
  titlee <- CMPS_hamby252_results$titlee[[i]]
  
  hamby252.plot.list[[1]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.max_scaled,
                       fill = as.factor(type_truth)), binwidth = 0.04) +
    labs(
      fill = "Comparison Type",
      x = expression(CMPS^"*"[max]),
      # title = expression(paste("Hamby252 - ", CMPS[max], " Distribution")),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = pretty(seq(0, 1, 0.05), n=10)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    font("x.text", size = 6) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  hamby252.plot.list[[2]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.maxbar_scaled,
                       fill = as.factor(type_truth)), binwidth = 0.04) +
    labs(
      x = expression(bar(CMPS^"*")[max]),
      fill = "Comparison Type",
      # title = expression(paste(
      #   "Hamby252 - ", bar(CMPS)[max], " Distribution"
      # )),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = pretty(seq(0, 1, 0.05), n=10)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  hamby252.plot.list[[3]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.diff_scaled,
                       fill = as.factor(type_truth)), binwidth = 0.05) +
    labs(
      x = expression(bar(CMPS^"*")[diff]),
      fill = "Comparison Type",
      # title = expression(paste(
      #   "Hamby252 - ", bar(CMPS)[max], " Distribution"
      # )),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = pretty(seq(0, 1, 0.05), n=10)) +
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

for (i in 1:length(CMPS_hamby252_results$signame)) {
  hamby252.cmps <- CMPS_hamby252_results$cmps.table[[i]]
  
  hamby252.plot.list <- list()
  titlee <- CMPS_hamby252_results$titlee[[i]]
  
  hamby252.plot.list[[1]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.max,
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
    geom_histogram(aes(x = cmps.maxbar,
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
    geom_histogram(aes(x = cmps.diff,
                       fill = as.factor(type_truth)), binwidth = 1) +
    labs(
      x = expression(bar(CMPS)[diff]),
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
CMPS_hamby252_results$plot[[4]]

CMPS_hamby252_results$plot_scaled[[1]]
CMPS_hamby252_results$plot_scaled[[4]]

##########
# variance ratio
table.cmps.diff <- lapply(CMPS_hamby252_results$cmps.table, function(tt) {
  r1 <- compute_var_ratio_anova(tt$cmps.diff, tt$type_truth, MS=FALSE)
  r2 <- compute_var_ratio_anova(tt$cmps.diff_scaled, tt$type_truth, MS=FALSE)
  r3 <- compute_var_ratio_anova(tt$cmps.max, tt$type_truth, MS=FALSE)
  r4 <- compute_var_ratio_anova(tt$cmps.max_scaled, tt$type_truth, MS=FALSE)
  r5 <- compute_var_ratio_anova(tt$cmps.maxbar, tt$type_truth, MS=FALSE)
  r6 <- compute_var_ratio_anova(tt$cmps.maxbar_scaled, tt$type_truth, MS=FALSE)
  c(r1, r2, r3, r4, r5, r6)
}) %>% do.call(rbind, .) %>% as.data.frame()

colnames(table.cmps.diff) <- c("cmps.diff", "cmps.diff_scaled",
                               "cmps.max", "cmps.max_scaled",
                               "cmps.maxbar", "cmps.maxbar_scaled")
table.cmps.diff$seg_length <- unlist(CMPS_hamby252_results$seg_length)

longer.tt <- table.cmps.diff %>% pivot_longer(cols = cmps.diff:cmps.maxbar_scaled,
                                 names_to = "score_name")
longer.tt <- longer.tt %>% mutate(
  scaled = str_detect(score_name, "scaled"),
  metric = str_remove(score_name, "_scaled")
)

longer.tt %>%
  ggplot(aes(x = seg_length, y = value)) +
  geom_line(aes(color = scaled, linetype = metric)) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
  theme_bw()

