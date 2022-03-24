## ----setup, echo=FALSE, message=FALSE, warning=FALSE----------------------
library(knitr)
library(kableExtra)
library(tidyverse)
library(patchwork)
library(ggpubr)

## ----func, echo=FALSE-----------------------------------------------------
bootstrap_k <- function(scores, k, K, value) {
  res <- replicate(K, {
    mean(sample(scores, size = k, replace = TRUE), na.rm = TRUE)
  })
  sum(res >= value) / K
}

compute_average_scores <- function(land1, land2, score, addNA = FALSE) {
  if (!is.numeric(land1)) land1 <- readr::parse_number(as.character(land1))
  if (!is.numeric(land2)) land2 <- readr::parse_number(as.character(land2))
  # assert_that(is.numeric(land1), is.numeric(land2), is.numeric(score))

  maxland <- max(land1, land2)
  fullframe <- data.frame(expand.grid(land1 = 1:maxland, land2 = 1:maxland))
  bcompare <- data.frame(land1, land2, score)

  fullframe <- fullframe %>% left_join(bcompare, by = c("land1", "land2"))

  fullframe <- fullframe %>% mutate(
    land1 = factor(land1, levels = 1:maxland),
    land2 = factor(land2, levels = 1:maxland)
  )
  # get averages, just in case
  matrix <- xtabs(score ~ land1 + land2,
    data = fullframe, addNA = addNA
  ) / xtabs(~land1 + land2, data = fullframe, addNA = addNA)

  matrix <- cbind(matrix, matrix)

  scores <- 1:maxland %>% sapply(FUN = function(i) {
    if (i == 1) {
      mean(diag(matrix), na.rm = TRUE)
    } else {
      i <- i - 1
      mean(diag(matrix[, -(1:i)]), na.rm = TRUE)
    }
  })
  scores
}

bullet_to_land_predict <- function(land1, land2, scores, difference, alpha = 0.05, addNA = FALSE) {
  if (!is.numeric(land1)) land1 <- readr::parse_number(as.character(land1))
  if (!is.numeric(land2)) land2 <- readr::parse_number(as.character(land2))
  avgs <- compute_average_scores(land2, land1, scores, addNA)

  p <- max(c(land1, land2))
  boot <- bootstrap_k(scores, p, 1000, max(avgs)) < alpha
  if (!is.numeric(boot)) {
    boot <- bootstrap_k(scores, p, 1000, max(avgs)) < alpha
    #     print("boot is not numeric")
    #      browser()
  }

  getdiff <- diff(sort(-avgs))[1]
  if (!is.numeric(getdiff)) print("getdiff is not numeric") # browser()
  if (getdiff > difference & boot) {
    # pick the maximum to determine the phase
    idx <- which.max(avgs)
    dd <- data.frame(
      land1,
      land2
    ) %>%
      mutate_if(function(.) !is.numeric(.), parse_number)
    dd$diff <- (dd$land1 - dd$land2) %% p + 1


    return(dd$diff == idx)
  } else {
    return(rep(FALSE, length = length(land1)))
  }
}


## ----bullet, echo=FALSE, out.width=".8\\textwidth", fig.cap="Photo of a traditionally rifled gun barrel (left) and a fired bullet (right)."----
# knitr::include_graphics("img/barrel_bullet_ps.png", dpi = 100)
knitr::include_graphics("ju-hofmann_files/figure-latex/barrel_bullet_ps.png", dpi = 144)


## ----process, echo=FALSE, out.width=".9\\textwidth", fig.cap="A framework of obtaining a bullet signature. (a) rendering from the 3D topographic scan of a land engraved area (LEA). The selected crosscut location is indicated by a thin white horizontal line. (b) view of the cross-section of the land engraved area at the white line in (a). (c) the crosscut data plotted in 2D; blue vertical lines indicate the position of left and right grooves. (d) the crosscut data after chopping the left and right grooves. (e) the fitted curvature using LOESS. (f) after removing the curvature from the crosscut data, the bullet signature is obtained"----
knitr::include_graphics("ju-hofmann_files/figure-latex/figure1_v2.png", dpi = 144)


## ---- eval=FALSE----------------------------------------------------------
## # install.packages("cmpsR")
## 
## library(cmpsR)
## data(bullets)
## 
## sig1 <- bullets$sigs[[2]]$sig
## sig2 <- bullets$sigs[[9]]$sig
## sig3 <- bullets$sigs[[10]]$sig
## 
## cmps.result.KM <- extract_feature_cmps(sig1, sig2)
## cmps.result.KNM <- extract_feature_cmps(sig1, sig3)


## ---- eval=FALSE----------------------------------------------------------
## extract_feature_cmps(
##   x,
##   y,
##   seg_length = 50,
##   Tx = 25,
##   npeaks.set = c(5, 3, 1),
##   include = NULL,
##   outlength = NULL
## )


## ---- eval=FALSE----------------------------------------------------------
## install.packages("cmpsR")


## ---- eval=FALSE----------------------------------------------------------
## # install.packages("remotes")
## remotes::install_github("willju-wangqian/cmpsR")


## ---- eval=TRUE-----------------------------------------------------------
library(cmpsR)
data(bullets)


## ----sigs, echo=FALSE, fig.cap="Signatures of all lands of bullet 1 in the top row, and of bullet 2 in the bottom row. Signatures in the second row are ordered to be in phase with the signatures above, i.e. matching signatures are displayed on top of each other."----
signatures <- bullets %>% unnest(sigs)
signatures <- signatures %>% mutate(
  bulletland = factor(bulletland, levels=c(paste(rep(c(1,2), each=6), c(1:6,2:6,1), sep="-")))
)
signatures %>% ggplot(aes(x = x/1000, y = sig)) + geom_line() + facet_wrap(~bulletland, ncol=6) +
  theme_bw() +
  xlab("Length in mm") +
  ylab("Relative height in micron")


## ---- eval=TRUE-----------------------------------------------------------
sigs1 <- bullets$sigs[bullets$bulletland == "2-5"][[1]]
sigs2 <- bullets$sigs[bullets$bulletland == "1-4"][[1]]

# compute cmps

# algorithm with multi-peak insepction at three different segment levels
cmps_with_multi_scale <- 
  extract_feature_cmps(sigs1$sig, sigs2$sig, 
                       npeaks.set = c(5,3,1), include = "full_result")

# algorithm with multi-peak inspection at the basis scale only
cmps_without_multi_scale <- 
  extract_feature_cmps(sigs1$sig, sigs2$sig, 
                       npeaks.set = 5, include = "full_result")


## ----alternative, echo = FALSE--------------------------------------------
lands <- unique(bullets$bulletland)

comparisons <- data.frame(expand.grid(land1 = lands[1:6], land2 = lands[7:12]), 
                          stringsAsFactors = FALSE)

comparisons <- comparisons %>% 
  left_join(bullets %>% select(bulletland, sig1=sigs),
            by = c("land1" = "bulletland")) %>%
  left_join(bullets %>% select(bulletland, sig2=sigs),
            by = c("land2" = "bulletland"))

comparisons <- comparisons %>% mutate(
  cmps = purrr::map2(sig1, sig2, .f = function(x, y) {
    extract_feature_cmps(x$sig, y$sig, include = "full")
    })
)

comparisons <- comparisons %>% 
  mutate(
    cmps_score = sapply(comparisons$cmps, function(x) x$CMPS.score),
    cmps_nseg = sapply(comparisons$cmps, function(x) x$nseg),
    cmps_scaled = cmps_score / cmps_nseg
  )

comparisons <- comparisons %>%
  mutate(
    bulletA = gsub("(\\d)-\\d", "\\1", land1),
    landA = gsub("\\d-(\\d)", "\\1", land1),
    bulletB = gsub("(\\d)-\\d", "\\1", land2),
    landB = gsub("\\d-(\\d)", "\\1", land2)
  )

dframe <- comparisons %>% select(-sig1, -sig2)

dframe$samesource <- with(dframe, bullet_to_land_predict(land1=landA, land2=landB, cmps_score, difference=1))


## ----sigplot, opts.label="codefig", echo=TRUE, out.width='\\textwidth', warning = FALSE----
sig.plot <- cmps_signature_plot(
  cmps_with_multi_scale
)
sig.plot$segment_shift_plot


## ----sigplot2, opts.label="codefig", echo=TRUE, out.width='\\textwidth', warning = FALSE----
sig.plot$signature_shift_plot


## ---- echo=TRUE, eval=TRUE------------------------------------------------
sig.plot$seg_shift


## ----segplot, opts.label="codefig", echo=TRUE, out.width='\\textwidth', warning = FALSE----
seg.plot <- cmps_segment_plot(
  cmps_with_multi_scale, 
  seg.idx = 6
)
seg.plot[[1]]$segment_plot


## ----segplot2, opts.label="codefig", echo=TRUE, out.width='\\textwidth', warning = FALSE----
seg.plot[[1]]$scale_ccf_plot


## ----segplot3, opts.label="codefig", echo=TRUE, out.width='\\textwidth', warning = FALSE----
library(ggpubr)

ggarrange(
  plotlist = 
    unlist(seg.plot, 
           recursive = FALSE),
  ncol = 2, 
  nrow = 3)


## ----tiles, echo=FALSE, fig.cap="CMPS scores of all 36 pairwise bullet signature comparisons for two bullets. Land engraving pairs generated by the same land (KM comparisons) are highlighted. Note that in this example the axis along Bullet 2 starts with Land 2. This corresponds to Phase 1 in equation (1).", out.width=".7\\textwidth", fig.align="center"----
dframe <- dframe %>% mutate(
  landA = paste0("L", landA),
  landB = paste0("L", landB),
  landB = factor(landB, levels = paste0("L", c(2:6,1))),
  bulletA = paste0("Bullet ", bulletA),
  bulletB = paste0("Bullet ", bulletB)
)

dframe %>% ggplot(aes(x = landA, y = landB, fill = cmps_score)) + 
  geom_tile() + 
  geom_tile(aes(colour="same land"), fill=NA, data = dframe %>% filter(samesource), size=1) + 
  scale_fill_gradient2("CMPS score", low = "gray80", high = "darkorange", midpoint = 6) + 
  scale_colour_manual("Source", values="darkorange") +
  facet_grid(bulletB ~ bulletA) + xlab("Bullet1 Lands") + 
  ylab("Bullet2 Lands") + 
  geom_text(aes(label=cmps_score)) +
  theme_bw() +
  theme(aspect.ratio = 1)


## ----tiles2, echo=FALSE, fig.cap="Plot (a) shows the highest possible CMPS scores (the total number of basis segments) for the 36 comparisons. (b) shows the scaled CMPS scores for the 36 comparisons.", out.width=".7\\textwidth", fig.align="center"----
p1 <- dframe %>% ggplot(aes(x = landA, y = landB)) + 
  geom_tile(fill = "grey") + 
  geom_tile(aes(colour="same land"), fill=NA, data = dframe %>% filter(samesource), size=1) + 
  scale_colour_manual("Source", values="darkorange") +
  facet_grid(bulletB ~ bulletA) + xlab("Bullet1 Lands") + 
  ylab("Bullet2 Lands") + 
  geom_text(aes(label=cmps_nseg)) +
  theme_bw() +
  theme(aspect.ratio = 1)

p2 <- dframe %>% ggplot(aes(x = landA, y = landB, fill = cmps_scaled)) + 
  geom_tile() + 
  geom_tile(aes(colour="same land"), fill=NA, data = dframe %>% filter(samesource), size=1) + 
  scale_fill_gradient2("Scaled CMPS", low = "gray80", high = "darkorange", midpoint = 0.3) + 
  scale_colour_manual("Source", values="darkorange") +
  facet_grid(bulletB ~ bulletA) + xlab("Bullet1 Lands") + 
  ylab("Bullet2 Lands") + 
  geom_text(aes(label=round(cmps_scaled, 2))) +
  theme_bw() +
  theme(aspect.ratio = 1)

ggarrange(p1, p2, 
          common.legend = TRUE, legend = "bottom", labels = c("(a)", "(b)"))


## ----load-results, eval=TRUE, echo=FALSE----------------------------------
all_results_container <- readRDS("data/CMPSpaper_results.rds")

CMPS_hamby252_results_seg <- all_results_container$h252_seg
CMPS_hamby252_results_npeak <- all_results_container$h252_npeak
CMPS_hamby44_results_seg <- all_results_container$h44_seg
CMPS_hamby44_results_npeak <- all_results_container$h44_npeak
rf.results <- all_results_container$rf

#### a global title for Hamby44 plots
com.title44 <- expression(paste(
  "Hamby 44 - ", CMPS[max], " and ", bar(CMPS)[max], " Distribution"
))

#### a global title for Hamby252 plots
com.title252 <- expression(paste(
  "Hamby 252 - ", CMPS[max], " and ", bar(CMPS)[max], " Distribution"
))


## ----rf-data-process, echo=FALSE------------------------------------------
rf.tt <- lapply(rf.results, function(rf_result) {
  with(rf_result, {
    tibble(rf.diff, rf.max, rf.maxbar) %>% 
      apply(2, compute_ss_ratio, label = type_truth, MS = FALSE)
  })
}) %>% do.call(rbind, .) %>% as_tibble()

colnames(rf.tt) <- c("cmps.diff", "cmps.max", "cmps.maxbar")

rf.tt$cmps.diff_scaled <- rf.tt$cmps.diff

rf.tt$study <- c("Hamby 252", "Hamby 44")
rf.tt <- rf.tt %>% pivot_longer(
  cols = c("cmps.diff", "cmps.max", "cmps.maxbar", "cmps.diff_scaled"),
  names_to = "Metric"
)

rf.tt$npeak <- "RF score* (Hare et al)"



## ----ss-original-cmps, echo=FALSE-----------------------------------------
result.container <- list()

result.container$npeak <- 
  c("5* (Chen et al)", "5-3-1* (Chen et al)", "5* (Chen et al)", "5-3-1* (Chen et al)")
result.container$Metric <- c("cmps.max", "cmps.max", "cmps.maxbar", "cmps.maxbar")
result.container$study <- "Hamby 252"

dd.list <- list()

dd.list[[1]] <- tibble(
  score = c(9:13, 16:24),
  freq = c(132, 282, 115, 19, 1, 1, 3,2,4,8,12,7,7,2),
  type = c(rep("KNM", 5), rep("KM", 9))
)

dd.list[[2]] <- tibble(
  score = c(1:5, 14:23),
  freq = c(7,343,168,24,7,1,1,2,3,6,7,8,8,7,3),
  type = c(rep("KNM", 5), rep("KM", 10))
)

dd.list[[3]] <- tibble(
  score = c(6:9, 13:19),
  freq = c(2,210,307,30,6,3,9,9,11,7,1),
  type = c(rep("KNM", 4), rep("KM", 7))
)

dd.list[[4]] <- tibble(
  score = c(0:2, 8:17),
  freq = c(36, 508, 5, 1,3,4,2,11,3,9,8,4,1),
  type = c(rep("KNM", 3), rep("KM", 10))
)

result.container$dd.list <- dd.list

result.container <- result.container %>% as_tibble() %>% 
  mutate(
    dd.table = purrr::map(dd.list, function(dd) {
      tibble(score = rep(dd$score, dd$freq),
             type_truth = rep(dd$type, dd$freq))
    }),
    value = purrr::map_dbl(dd.table, function(dd) {
      compute_ss_ratio(dd$score, dd$type_truth, MS = FALSE)
    })
  ) 


## ---- eval=FALSE----------------------------------------------------------
## extract_feature_cmps(
##   x, y,
##   seg_length = 50,
##   Tx = 25,
##   npeaks.set = c(5,3,1),
##   include = "nseg"
## )


## ---- echo=FALSE, eval=FALSE----------------------------------------------
## with(CMPS_hamby252_results_npeak$cmps.table[[3]], {
##   tibble(cmps.max, cmps.maxbar) %>%
##     apply(2, compute_ss_ratio, label = type_truth, MS = FALSE)
## })
## 
## with(CMPS_hamby252_results_seg$cmps.table[[2]], {
##   tibble(cmps.max, cmps.maxbar, cmps.diff, cmps.diff_scaled) %>%
##     apply(2, compute_ss_ratio, label = type_truth, MS = FALSE)
## })
## 


## ----color-settings, echo=FALSE-------------------------------------------

cols_1 <- c("#748BA7", "#C583AE","#FFDFAA", "#D9F0A0")
cols_2 <- c("#4C688B", "#A45287","#D4AC6A", "#ABC864")
# #####  Color Palette by Paletton.com
# #####  Palette URL: http://paletton.com/#uid=7561m0kllllaFw0g0qFqFg0w0aF
# 
# 
# *** Primary color:
# 
#    shade 0 = #832C65 = rgb(131, 44,101) = rgba(131, 44,101,1) = rgb0(0.514,0.173,0.396)
#    shade 1 = #C583AE = rgb(197,131,174) = rgba(197,131,174,1) = rgb0(0.773,0.514,0.682)
#    shade 2 = #A45287 = rgb(164, 82,135) = rgba(164, 82,135,1) = rgb0(0.643,0.322,0.529)
#    shade 3 = #621046 = rgb( 98, 16, 70) = rgba( 98, 16, 70,1) = rgb0(0.384,0.063,0.275)
#    shade 4 = #42002B = rgb( 66,  0, 43) = rgba( 66,  0, 43,1) = rgb0(0.259,0,0.169)
# 
# *** Secondary color (1):
# 
#    shade 0 = #AA7F39 = rgb(170,127, 57) = rgba(170,127, 57,1) = rgb0(0.667,0.498,0.224)
#    shade 1 = #FFDFAA = rgb(255,223,170) = rgba(255,223,170,1) = rgb0(1,0.875,0.667)
#    shade 2 = #D4AC6A = rgb(212,172,106) = rgba(212,172,106,1) = rgb0(0.831,0.675,0.416)
#    shade 3 = #805715 = rgb(128, 87, 21) = rgba(128, 87, 21,1) = rgb0(0.502,0.341,0.082)
#    shade 4 = #553500 = rgb( 85, 53,  0) = rgba( 85, 53,  0,1) = rgb0(0.333,0.208,0)
# 
# *** Secondary color (2):
# 
#    shade 0 = #2B4970 = rgb( 43, 73,112) = rgba( 43, 73,112,1) = rgb0(0.169,0.286,0.439)
#    shade 1 = #748BA7 = rgb(116,139,167) = rgba(116,139,167,1) = rgb0(0.455,0.545,0.655)
#    shade 2 = #4C688B = rgb( 76,104,139) = rgba( 76,104,139,1) = rgb0(0.298,0.408,0.545)
#    shade 3 = #143054 = rgb( 20, 48, 84) = rgba( 20, 48, 84,1) = rgb0(0.078,0.188,0.329)
#    shade 4 = #051B38 = rgb(  5, 27, 56) = rgba(  5, 27, 56,1) = rgb0(0.02,0.106,0.22)
# 
# *** Complement color:
# 
#    shade 0 = #81A035 = rgb(129,160, 53) = rgba(129,160, 53,1) = rgb0(0.506,0.627,0.208)
#    shade 1 = #D9F0A0 = rgb(217,240,160) = rgba(217,240,160,1) = rgb0(0.851,0.941,0.627)
#    shade 2 = #ABC864 = rgb(171,200,100) = rgba(171,200,100,1) = rgb0(0.671,0.784,0.392)
#    shade 3 = #5B7814 = rgb( 91,120, 20) = rgba( 91,120, 20,1) = rgb0(0.357,0.471,0.078)
#    shade 4 = #395000 = rgb( 57, 80,  0) = rgba( 57, 80,  0,1) = rgb0(0.224,0.314,0)
# 
# 
# #####  Generated by Paletton.com (c) 2002-2014


## ----plot_computation, echo=FALSE, fig.width=8, fig.height = 8, out.width='\\linewidth'----
# 252 - seg
table.252.cmps_seg <- lapply(CMPS_hamby252_results_seg$cmps.table, function(tt) {
  type_truth <- tt$type_truth
  tt %>% select(cmps.diff:cmps.maxbar_scaled) %>% 
    apply(2, compute_ss_ratio, label=type_truth, MS=FALSE)
}) %>% do.call(rbind, .) %>% as.data.frame()

table.252.cmps_seg$seg_length <- unlist(CMPS_hamby252_results_seg$seg_length)

table.252.cmps_seg.long <- table.252.cmps_seg %>%
  pivot_longer(cols = c(cmps.diff, cmps.diff_scaled, 
                        cmps.max, cmps.maxbar), 
               names_to = "Metric")

table.252.cmps_seg.long <- table.252.cmps_seg.long %>% 
  mutate(study = "Hamby 252")


labels <- expression( CMPS[max], bar(CMPS)[diff], {bar(CMPS)^symbol("*")}[diff], bar(CMPS)[max])

p.seg_252_plot <- table.252.cmps_seg.long %>%
  ggplot(aes(x = seg_length +3*as.numeric(factor(Metric))-9)) +
  geom_segment(aes(xend =  seg_length +3*as.numeric(factor(Metric))-9,
               y = value, yend = 0, colour = Metric)) + 
  geom_point(aes(y =  value, color = Metric), 
             alpha = 0.8,
             position=position_dodge(width = 0.5), size = 3.5) +
  scale_x_continuous(breaks = table.252.cmps_seg$seg_length, 
                     minor_breaks = NULL) +
  theme_bw() +
  # theme(legend.position="none") +
  ggtitle("Results for Hamby 252") +
  ylab("Sum of Squares Ratio") + 
  xlab("Segment Length") +
  scale_colour_manual(values=cols_2, labels = labels) +
  scale_fill_manual(values=cols_1, labels = labels) 

# 44 - seg
table.44.cmps_seg <- lapply(CMPS_hamby44_results_seg$cmps.table, function(tt) {
  type_truth <- tt$type_truth
  tt %>% select(cmps.diff:cmps.maxbar_scaled) %>% 
    apply(2, compute_ss_ratio, label=type_truth, MS=FALSE)
}) %>% do.call(rbind, .) %>% as.data.frame()

table.44.cmps_seg$seg_length <- unlist(CMPS_hamby44_results_seg$seg_length)

table.44.cmps_seg.long <- table.44.cmps_seg %>%
  pivot_longer(cols = c(cmps.diff, cmps.diff_scaled, 
                        cmps.max, cmps.maxbar),
               names_to = "Metric")

table.44.cmps_seg.long <- table.44.cmps_seg.long %>% mutate(
  study = "Hamby 44"
)

table.cmps_seg.long <- bind_rows(table.44.cmps_seg.long,
                                 table.252.cmps_seg.long)

breaks_fun <- function(x) {
  if(x[2] < 220) {
    c(25,  50,  75, 100, 125, 150, 175, 200)
  } else {
    c(30, 45,  61,  90, 122, 150, 180, 210)
  }
}

p.seg_plot <- table.cmps_seg.long %>%
  mutate(
    Metric = factor(Metric,
                    levels = c("cmps.max", "cmps.diff", "cmps.diff_scaled", "cmps.maxbar"))
  ) %>% 
  ggplot(aes(x = seg_length+3*as.numeric((Metric))-9)) +
  geom_segment(aes(xend = seg_length+3*as.numeric((Metric))-9, 
                   y = value, yend= 0, color=Metric),
               position=position_dodge(width = 0.6)) +
  geom_point(aes(y =  value, color = Metric), 
             alpha = 0.8,
             position=position_dodge(width = 0.5), size = 3.5) +
  facet_grid(.~study, scales = "free_x") +
  scale_x_continuous(breaks = breaks_fun, 
                     minor_breaks = NULL) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylab("Sums of Squares Ratio") + 
  xlab("Segment length") +
  scale_colour_manual(values=cols_2, labels = labels) +
  scale_fill_manual(values=cols_1, labels = labels) 

# npeak - 252
table.252.cmps_npeak <- lapply(CMPS_hamby252_results_npeak$cmps.table, function(tt) {
  type_truth <- tt$type_truth
  tt %>% select(cmps.diff:cmps.maxbar_scaled) %>% 
    apply(2, compute_ss_ratio, label=type_truth, MS=FALSE)
}) %>% do.call(rbind, .) %>% as.data.frame()

table.252.cmps_npeak$npeak <- CMPS_hamby252_results_npeak$npeaks.set %>% 
  purrr::map(.f = paste, collapse = "-") %>% sapply(as.character)

table.252.cmps_npeak.longer <- table.252.cmps_npeak %>% as_tibble() %>% 
  filter(npeak %in% c("5", "5-3-1", "10-8-6-4-2-1", "10-7-4-2")) %>% 
  mutate(
    npeak = factor(
      npeak,
      levels = npeak[order(cmps.diff_scaled, decreasing = FALSE)])
  ) %>% 
  pivot_longer(cols = c(cmps.diff, cmps.diff_scaled, cmps.max, cmps.maxbar),
               names_to = "Metric") 
table.252.cmps_npeak.longer$study <- "Hamby 252"

# npeak - 44
table.44.cmps_npeak <- lapply(CMPS_hamby44_results_npeak$cmps.table, function(tt) {
  type_truth <- tt$type_truth
  tt %>% select(cmps.diff:cmps.maxbar_scaled) %>% 
    apply(2, compute_ss_ratio, label=type_truth, MS=FALSE)
}) %>% do.call(rbind, .) %>% as.data.frame()

table.44.cmps_npeak$npeak <- CMPS_hamby44_results_npeak$npeaks.set %>% 
  purrr::map(.f = paste, collapse = "-") %>% sapply(as.character)

table.44.cmps_npeak.longer <- table.44.cmps_npeak %>% as_tibble() %>% 
  filter(npeak %in% c("5", "5-3-1", "10-8-6-4-2-1", "10-7-4-2")) %>% 
  mutate(
    npeak = factor(
      npeak,
      levels = npeak[order(cmps.diff_scaled, decreasing = FALSE)])
  ) %>% 
  pivot_longer(cols = c(cmps.diff, cmps.diff_scaled, cmps.max, cmps.maxbar),
               names_to = "Metric") 
table.44.cmps_npeak.longer$study <- "Hamby 44"

selected_var <- c("npeak", "Metric", "study", "value")

table.cmps_neapk.long <- bind_rows(
  table.252.cmps_npeak.longer %>% select(all_of(selected_var)),
  table.44.cmps_npeak.longer %>% select(all_of(selected_var)),
  rf.tt %>% select(all_of(selected_var)),
  result.container %>% select(all_of(selected_var))
)

order.level <- table.cmps_neapk.long %>% 
  filter(Metric == "cmps.maxbar", study == "Hamby 252") %>%
  {.$npeak[order(.$value, decreasing = FALSE)]}

table.cmps_neapk.long <- table.cmps_neapk.long %>% mutate(
  npeak = factor(npeak, levels = order.level)
)

npeak_convert <- with(table.cmps_neapk.long, {
  tt <- tibble(npeak = levels(npeak))
  tt$xstart <- seq.int(from = 0, by = 20, length.out = length(tt$npeak))
  tt
})

table.cmps_neapk.long <- table.cmps_neapk.long %>% left_join(
  npeak_convert,
  by = "npeak"
)

table.cmps_neapk.long <- table.cmps_neapk.long %>% mutate(
  rf_status = if_else(npeak == "RF score* (Hare et al)",
                      TRUE, FALSE)
)

table.cmps_neapk.long$Metric <- 
  table.cmps_neapk.long$Metric %>% 
  factor(levels = c("cmps.max", "cmps.diff", "cmps.diff_scaled", "cmps.maxbar"))

p.npeak_plot <- table.cmps_neapk.long %>%
  ggplot(aes(x = xstart + 3*as.numeric(factor(Metric))-9)) +
  geom_segment(aes(xend = xstart + 3*as.numeric(factor(Metric))-9, 
                   y = value, yend= 0, color=Metric),
               position=position_dodge(width = 0.6)) +
  geom_point(aes(y =  value, color = Metric, shape = rf_status), 
             alpha = 0.8,
             position=position_dodge(width = 0.5), size = 3.5) +
  facet_grid(.~study) +
  ylab("Sums of Squares Ratio") + 
  xlab("Number of Peaks") +
  coord_flip() +
  scale_x_continuous(breaks = npeak_convert$xstart, labels = npeak_convert$npeak ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(shape = "none") + 
  scale_shape_manual(values = c(16, 1)) +
  scale_colour_manual(values=cols_2, labels = labels) +
  scale_fill_manual(values=cols_1, labels = labels) 


## ----result1_252, echo=FALSE, out.width="400px", fig.cap="Distribution of $\\mathrm{CMPS_{max}}$ and $\\mathrm{\\overline{CMPS}_{max}}$ for Hamby 252; outliers are removed in bullet signatures; \\texttt{seg\\_length = 50}, \\texttt{Tx = 25}, \\texttt{npeaks.set = c(5,3,1)} ", fig.height=3----
# [4] "span1=0.25, span2=0.03, rm_outliers = TRUE, 5025, 531"  
# knitr::include_graphics("img/hamby252_v1.png", dpi = 100)
# CMPS_hamby252_results$plot[[2]]

p1 <- metric_plot_helper(CMPS_hamby252_results_seg$cmps.table[[2]], "cmps.max", 
                         subtitle = "npeaks.set: 5-3-1, seg_length: 50") +
  xlab(labels[1])
p2 <- metric_plot_helper(CMPS_hamby252_results_seg$cmps.table[[2]], "cmps.maxbar",
                         subtitle = "npeaks.set: 5-3-1, seg_length: 50") +
  xlab(labels[4])

p.result1_252 <- ggarrange(plotlist = list(p1, p2), nrow = 1, ncol = 2,
               common.legend = TRUE, legend = "bottom")
p.result1_252 <- annotate_figure(p.result1_252, top = text_grob(com.title252))
p.result1_252


## ---- eval=FALSE----------------------------------------------------------
## extract_feature_cmps(
##   x, y,
##   seg_length = 61,
##   Tx = 30,
##   npeaks.set = c(5,3,1),
##   include = "nseg"
## )


## ----result1_44, echo=FALSE, out.width="400px", fig.cap="Distribution of $\\mathrm{CMPS_{max}}$ and $\\mathrm{\\overline{CMPS}_{max}}$ for Hamby 44; outliers are removed in bullet signatures; \\texttt{seg\\_length = 61}, \\texttt{Tx = 30}, \\texttt{npeaks.set = c(5,3,1)} ", fig.height=3, warning=FALSE----
p1 <- metric_plot_helper(CMPS_hamby44_results_seg$cmps.table[[3]], "cmps.max",
                         subtitle = "npeaks.set: 5-3-1, seg_length: 61") +
  xlab(labels[1])
p2 <- metric_plot_helper(CMPS_hamby44_results_seg$cmps.table[[3]], "cmps.maxbar",
                         subtitle = "npeaks.set: 5-3-1, seg_length: 61") +
  xlab(labels[4])
p.result1_44 <- ggarrange(plotlist = list(p1, p2), nrow = 1, ncol = 2,
               common.legend = TRUE, legend = "bottom")
p.result1_44 <- annotate_figure(p.result1_44, top = text_grob(com.title44))
p.result1_44


## ----param_seg_plot, echo=FALSE, out.width="400px", fig.cap="Comparison of results from the CMPS algorithm based on different basis segment lengths (parameter \\texttt{seg\\_length}). Only the $\\mathrm{CMPS_{max}}$ metric suggests that the default values for the basis segment length result in the best separation. Better separation is achieved based on the modified CMPS metrics, inclusing the newly suggested ones. For Hamby 252 these metrics agree on a segment length of 75, and a segment length of 122 for Hamby 44 yields better results", fig.height=3, warning = FALSE----
p.seg_plot 


## ----param_npeak_plot, echo=FALSE, out.width="400px", fig.cap="Comparison of CMPS results based on different strategies of number of peak selections. Starred results compare CMPS performance with results published in the literature. Results for the random forest score are represented with circles because the metrics are computed not based on the CMPS scores, but the random forest scores with the same logic. Since random forest scores lie within the interval $[0, 1]$, scaling the random forest scores will not change the results.", fig.height=3, warning = FALSE----
p.npeak_plot
