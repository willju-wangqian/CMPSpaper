library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(CMPS)
library(ggpubr)
library(parallel)
library(assertthat)

source("code/func_collection.R")

#########################
# run the example in the chapter
data("bullets")
bullets$bulletland

lands <- unique(bullets$bulletland)
comparisons <- data.frame(expand.grid(land1 = lands, land2 = lands), stringsAsFactors = FALSE)

comparisons <- comparisons %>% mutate(aligned = purrr::map2(.x = land1, .y = land2, 
                                                            .f = function(xx, yy) {
                                                              land1 <- bullets$sigs[bullets$bulletland == xx][[1]]
                                                              land2 <- bullets$sigs[bullets$bulletland == yy][[1]]
                                                              land1$bullet <- "first-land"
                                                              land2$bullet <- "second-land"
                                                              
                                                              # bulletxtrctr::sig_align(land1$sig, land2$sig)
                                                              sig_align(land1$sig, land2$sig)
                                                            }))

purrr::map_dbl(comparisons$aligned, .f = function(x) x$lag)

subset(comparisons, land1 == "2-3" & land2 == "1-2")$aligned[[1]]$lands %>% 
  mutate(`b2-l3` = sig1, `b1-l2` = sig2) %>% select(-sig1, -sig2) %>% tidyr::gather(sigs, 
                                                                                    value, `b2-l3`, `b1-l2`) %>% ggplot(aes(x = x, y = value, colour = sigs)) + 
  geom_line() + theme_bw() + scale_color_brewer(palette = "Dark2")

comparisons <- comparisons %>% mutate(ccf0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_ccf(x$lands)), 
                                      lag0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_lag(x$lands)), 
                                      D0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_D(x$lands)), 
                                      length0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_length(x$lands)), 
                                      overlap0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_overlap(x$lands)), 
                                      striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75), 
                                      cms_per_mm = purrr::map2(striae, 
                                                               aligned, .f = function(s, a) {
                                                                 extract_feature_cms_per_mm(s$lines, a$lands, resolution = 1.5625)
                                                               }), 
                                      matches0 = striae %>% purrr::map_dbl(.f = function(s) {
                                        bulletxtrctr:::extract_helper_feature_n_striae(s$lines, type = "peak", 
                                                                                       match = TRUE)
                                      }), 
                                      mismatches0 = striae %>% purrr::map_dbl(.f = function(s) {
                                        bulletxtrctr:::extract_helper_feature_n_striae(s$lines, type = "peak", 
                                                                                       match = FALSE)
                                      }), 
                                      bulletA = gsub("([1-2])-([1-6])", "\\1", land1), 
                                      bulletB = gsub("([1-2])-([1-6])", 
                                                     "\\1", land2), 
                                      landA = gsub("([1-2])-([1-6])", "\\2", land1), 
                                      landB = gsub("([1-2])-([1-6])", 
                                                   "\\2", land2))

########################
# test on random forest code
tmp.lands <-
  c(
    b252 %>% filter(bullet_id %in% b.cb[, cb.idx][1]) %>% pull(scan_id),
    b252 %>% filter(bullet_id %in% b.cb[, cb.idx][2]) %>% pull(scan_id)
  )
tmp.comp <-
  data.frame(expand.grid(land1 = tmp.lands[1:6], land2 = tmp.lands[7:12]),
             stringsAsFactors = FALSE)
tmp.comp$landidx1 <- tmp.comp$land1 %>% as.character() %>% str_sub(-1, -1) %>% 
  as.numeric()
tmp.comp$landidx2 <- tmp.comp$land2 %>% as.character() %>% str_sub(-1, -1) %>% 
  as.numeric()

tmp.comp <- tmp.comp %>% 
  mutate(
    aligned = purrr::map2(.x = land1, .y = land2, 
                          .f = function(xx, yy) {
                            land1 <- b252$sigs_main[b252$scan_id == xx][[1]]
                            land2 <- b252$sigs_main[b252$scan_id == yy][[1]]
                            land1$bullet <- "first-land"
                            land2$bullet <- "second-land"
                            
                            sig_align(land1$sig, land2$sig)
                            # bulletxtrctr::sig_align(land1$sig, land2$sig)
                          }))

tmp.comp$aligned %>% purrr::map_dbl(.f = function(x) x$lag)
tmp.comp$aligned %>% purrr::map_dbl(.f = function(x) my_extract_feature_lag(x$lands))
tmp.comp$aligned %>% purrr::map_dbl(.f = function(x) bulletxtrctr::extract_feature_lag(x$lands))
tmp.comp$aligned[[1]]$lands %>% extract_feature_lag()
tmp.comp$aligned[[2]]$lands %>% extract_feature_lag()


tmp.comp <- tmp.comp %>% 
  mutate(
    ccf0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_ccf(x$lands)), 
    lag0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_lag(x$lands)), 
    D0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_D(x$lands)), 
    length0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_length(x$lands)), 
    overlap0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_overlap(x$lands)), 
    striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75), 
    cms_per_mm = purrr::map2(striae, aligned, .f = function(s, a) {
      extract_feature_cms_per_mm(s$lines, a$lands, resolution = 1.5625)
    }), 
    matches0 = striae %>% purrr::map_dbl(.f = function(s) {
      bulletxtrctr:::extract_helper_feature_n_striae(s$lines, type = "peak", 
                                                     match = TRUE)
    }), 
    mismatches0 = striae %>% purrr::map_dbl(.f = function(s) {
      bulletxtrctr:::extract_helper_feature_n_striae(s$lines, type = "peak", 
                                                     match = FALSE)
    }), 
    bulletA = gsub("([1-2])-([1-6])", "\\1", land1), 
    bulletB = gsub("([1-2])-([1-6])", "\\1", land2), 
    landA = gsub("([1-2])-([1-6])", "\\2", land1), 
    landB = gsub("([1-2])-([1-6])", "\\2", land2))

tmp.comp <- tmp.comp %>% 
  mutate(
    striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75),
    features = purrr::map2(.x = aligned, .y = striae, 
                           .f = my_extract_feature_all, resolution = 1.5625), 
    legacy_features = purrr::map(striae, 
                                 extract_features_all_legacy, resolution = 1.5625)) %>% 
  tidyr::unnest(legacy_features)

require(randomForest)

tmp.comp$rfscore <- predict(bulletxtrctr::rtrees, newdata = tmp.comp, 
                               type = "prob")[, 2]

with(tmp.comp, {
  compute_avgscore_denoise(landidx1, landidx2, rfscore)
})

tmp.comp %>% 
  ggplot(aes(x = land1, y = land2, fill = rfscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "grey80", high = "darkorange", 
                       midpoint = .5) +
  # facet_grid(bulletB~bulletA, labeller = "label_both") +
  xlab("Land A") +
  ylab("Land B") +
  theme(aspect.ratio = 1)

tmp.comp$aligned[[1]]$lands %>% extract_feature_lag()

tmp.comp$aligned[[1]]$lands

library(CMPS)
data("bullets")
sig1 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]$sig
sig2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]$sig

segments <- get_segs(sig1, len = 50)
seg1 <- segments$segs[[1]]

bulletxtrctr::get_ccf(seg1, sig1, min.overlap = 50) -> r1
CMPS::get_ccf4(sig1, seg1, min.overlap = 50) -> r2

identical(r1$lag, r2$lag)
r1$lag %>% head()
r2$lag %>% tail()

sig_align(sig1, sig2) -> r1
bulletxtrctr::sig_align(sig1, sig2) -> r2
identical(r1$lands$sig1, r2$lands$sig1)
identical(r1$lands, r2$lands)
identical(r1$ccf)
r1$lands %>% head()
r2$lands %>% head()
