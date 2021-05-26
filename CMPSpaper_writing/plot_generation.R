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
# check 1-4 vs 2-5, segment 8




