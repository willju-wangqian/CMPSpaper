library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(CMPS)
library(ggpubr)
library(parallel)

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
                           .f = extract_features_all, resolution = 1.5625), 
    legacy_features = purrr::map(striae, 
                                 extract_features_all_legacy, resolution = 1.5625)) %>% 
  tidyr::unnest(legacy_features)


tmp.comp$aligned[[1]]$lands %>% extract_feature_lag()

tmp.comp$aligned[[1]]$lands



sig1 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]$sig
sig2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]$sig

get_ccf2(sig1, sig2) -> r1
get_ccf4(sig1, sig2) -> r2

sig_align(sig1, sig2) -> r1
bulletxtrctr::sig_align(sig1, sig2) -> r2
identical(r1$lands$sig1, r2$lands$sig1)
identical(r1$lands, r2$lands)
identical(r1$ccf)
r1$lands %>% head()
r2$lands %>% head()
