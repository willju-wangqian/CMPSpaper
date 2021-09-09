library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(CMPS)
library(ggpubr)
library(parallel)

source("code/func_collection.R")

###########
# from the x3p data, re-run CMPS results for hamby 44
# hamby44_results <-
#   readRDS("C:/Research/Bullet/Bullet_Rcode/hamby44_results.rds")

# b44.full <-
#   readRDS("C:/Research/Bullet/Bullet_Rcode/hamby44/hamby44_manual5.rds")
b44.full <- 
  readRDS("~/Research/CMPSpaper/preconsideration/hamby44_manual5.rds")

# rds file that stores indices of damaged scans
# idf.idx <-
#   readRDS("C:/Research/Bullet/Bullet_Rcode/hamby44/idfidx.rds")
# idf.idx <- 
#   readRDS("~/Research/CMPSpaper/preconsideration/idfidx.rds")
# 
# idf.scan_id <- b44.full %>% slice(idf.idx) %>% .$scan_id
idf.idx <- c(23,35,39,41,102,126,149)
idf.scan_id <- b44.full %>% slice(idf.idx) %>% .$scan_id

# generate bullet_id
bulletid.tb <- b44.full %>% select(scan_id)
bulletid.tb <- bulletid.tb %>% mutate(
  barrel_id = sapply(strsplit(scan_id, "-"), "[[", 1),
  bullet = sapply(strsplit(scan_id, "-"), "[[", 2),
  land_id = sapply(strsplit(scan_id, "-"), "[[", 3),
  bullet_id = paste(barrel_id, bullet)
)


# remove x3p files
b44 <- b44.full %>% select(-x3p,-crosscut,-grooves,-ccdata)

# get all comparisons
b44$bullet_id <- bulletid.tb %>% pull(bullet_id)
b.cb <- unique(b44$bullet_id) %>% utils::combn(m = 2)
p <- 1:595

# get ground-truth
key.table <- tibble(
  bullet.id = c(
    "K",
    "O",
    "L",
    "P",
    "J",
    "H",
    "I",
    "Y",
    "G",
    "E",
    "U",
    "X",
    "F",
    "T",
    "S"
  ),
  barrel.id = c(1, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 9, 10)
)

tmp.table.44 <- tibble(bullet1 = b.cb[1, ],
                       bullet2 = b.cb[2, ])

tmp.table.44$b1t <-
  tmp.table.44$bullet1 %>% str_extract("'[A-Z]'") %>% str_extract("[A-Z]")
tmp.table.44$b2t <-
  tmp.table.44$bullet2 %>% str_extract("'[A-Z]'") %>% str_extract("[A-Z]")

tmp.table.44$b1b <-
  tmp.table.44 %>% left_join(key.table, by = c("b1t" = "bullet.id")) %>% pull(barrel.id)
tmp.table.44$b2b <-
  tmp.table.44 %>% left_join(key.table, by = c("b2t" = "bullet.id")) %>% pull(barrel.id)

tmp.table.44 <- tmp.table.44 %>% mutate(
  b1b = if_else(
    is.na(b1b),
    bullet1 %>% str_extract("Barrel \\d+") %>% str_extract("\\d+") %>% as.numeric(),
    b1b
  ),
  b2b = if_else(
    is.na(b2b),
    bullet2 %>% str_extract("Barrel \\d+") %>% str_extract("\\d+") %>% as.numeric(),
    b2b
  ),
  type_truth = if_else(b1b == b2b, "KM", "KNM")
)

tmp.table.44 <- tmp.table.44 %>% mutate(
  type = if_else(
    bullet1 %>% str_extract("(Unknowns|Barrel) (\\d+|[A-Z])") ==
      bullet2 %>% str_extract("(Unknowns|Barrel) (\\d+|[A-Z])"),
    "KM",
    "KNM"
  ),
  type = if_else(
    bullet1 %>% str_extract("(Unknowns|Barrel) (\\d+|[A-Z])") == "Unknowns B" |
      bullet2 %>% str_extract("(Unknowns|Barrel) (\\d+|[A-Z])") == "Unknowns B",
    "UKN",
    type
  )
)


N <- 9
CMPS_hamby44_results <- list()
CMPS_hamby44_results$span1 <- as.list(rep(0.25, N))
CMPS_hamby44_results$signame <-
  list("n5", "n51", "n52", "n531", "n10-62", "n6421",
       "n10-742", "n10-6421", "n10-86421")
  # list("s30", "s45", "s61", "s90", "s122", "s150", "s180", "s210")
CMPS_hamby44_results$npeaks.set <-
  list(c(5),
       c(5,1),
       c(5,2),
       c(5,3,1),
       c(10,6,2),
       c(6,4,2,1),
       c(10,7,4,2),
       c(10,6,4,2,1),
       c(10, 8, 6, 4, 2, 1))
  # lapply(1:N, function(t) c(5,3,1))
CMPS_hamby44_results$seg_length <- as.list(rep(61, N))
  # as.list(c(30, 45, 61, 90, 122, 150, 180, 210))
CMPS_hamby44_results$outlength <- vector(mode = "list", length = N)
CMPS_hamby44_results$outlength[[N]] <- c(61, 86, 111, 136, 161, 186)
CMPS_hamby44_results$Tx <- as.list(rep(30, N))
CMPS_hamby44_results$titlee <- list()
CMPS_hamby44_results$filename <- list()
for (i in 1:N) {
  CMPS_hamby44_results$titlee[[i]] <-
    paste0(
      "npeaks.set=c(",
      paste(CMPS_hamby44_results$npeaks.set[[i]], collapse = ","),
      ")",
      ", len=",
      CMPS_hamby44_results$seg_length[[i]],
      ", Tx=",
      CMPS_hamby44_results$Tx[[i]],
      ", \nspan1=",
      CMPS_hamby44_results$span1[[i]],
      ", span2=0.03"
    )
  CMPS_hamby44_results$filename[[i]] <- 
    paste(
      "hamby44",
      CMPS_hamby44_results$span1[[i]]*100,
      paste(CMPS_hamby44_results$npeaks.set[[i]], collapse = "-"),
      CMPS_hamby44_results$seg_length[[i]],
      CMPS_hamby44_results$Tx[[i]],
      sep = "_"
    )
}
CMPS_hamby44_results$cmps.table <- list()
CMPS_hamby44_results$plot <- list()

com.title44 <- expression(paste(
  "Hamby 44 - ", CMPS[max], " and ", bar(CMPS)[max], " Distribution"
))


# start of the loop!
for (i in 1:N) {
  # i <- 1
  if((i == 1) || (CMPS_hamby44_results$span1[[i]] != CMPS_hamby44_results$span1[[i-1]])){
    b44[[CMPS_hamby44_results$signame[[i]]]] <- purrr::map2(
      .x = b44.full$ccdata,
      .y = b44.full$grooves,
      .f = function(x, y) {
        cc_get_signature(
          ccdata = x,
          grooves = y,
          span1 = CMPS_hamby44_results$span1[[i]],
          span2 = 0.03
        )
      })
  } else {
    b44[[CMPS_hamby44_results$signame[[i]]]] <- b44[[CMPS_hamby44_results$signame[[i-1]]]]
  }
}
  
for (i in 1:N) {
  # process of removing outliers
  tt <-
    lapply(b44[[CMPS_hamby44_results$signame[[i]]]], function(x) {
      x$sig
    }) %>% unlist()
  qtt <- quantile(tt, na.rm = TRUE)
  multi.iqr <- 3
  outrange <-
    c(qtt[2] - multi.iqr * (qtt[4] - qtt[2]),
      qtt[4] + multi.iqr * (qtt[4] - qtt[2]))
  
  b44 <- b44 %>% mutate(sigs_main = purrr::map(
    .x = .data[[CMPS_hamby44_results$signame[[i]]]],
    .f = function(x) {
      tmp <- x$sig
      tmp[tmp < outrange[1] | tmp > outrange[2]] <- NA
      list(sig = tmp)
    }
  ))
  
  ## parallel setup
  # cl <- makeCluster(7)
  # par.setup <- parLapply(cl, 1:length(cl),
  #                        function(xx) {
  #                          library(tidyverse)
  #                          library(bulletxtrctr)
  #                          library(x3ptools)
  #                          library(CMPS)
  #                        })
  # clusterExport(cl,
  #               c("b44", "b.cb", "p",
  #                 "CMPS_hamby44_results", "i"),
  #               envir = globalenv())
  
  
  # start the parallel computing
  system.time({
    # tmp.44.list <- parLapply(cl, p, function(cb.idx) {
    tmp.44.list <- mclapply(p, function(cb.idx) {
      
      tmp.lands <-
        c(
          b44 %>% filter(bullet_id %in% b.cb[, cb.idx][1]) %>% pull(scan_id),
          b44 %>% filter(bullet_id %in% b.cb[, cb.idx][2]) %>% pull(scan_id)
        )
      tmp.comp <-
        data.frame(expand.grid(land1 = tmp.lands[1:6], land2 = tmp.lands[7:12]),
                   stringsAsFactors = FALSE)
      tmp.comp$landidx1 <- tmp.comp$land1 %>% as.character() %>% str_sub(-1, -1) %>% 
        as.numeric()
      tmp.comp$landidx2 <- tmp.comp$land2 %>% as.character() %>% str_sub(-1, -1) %>% 
        as.numeric()
      
      
      # cat(cb.idx, "- 1; ")
      tmp.comp$cmps <-
        purrr::map2(
          .x = tmp.comp$land1,
          .y = tmp.comp$land2,
          .f = function(xx, yy) {
            sig1 <- b44$sigs_main[b44$scan_id == xx][[1]]$sig
            sig2 <- b44$sigs_main[b44$scan_id == yy][[1]]$sig
            
            extract_feature_cmps(
              sig1,
              sig2,
              npeaks.set = CMPS_hamby44_results$npeaks.set[[i]],
              seg_length = CMPS_hamby44_results$seg_length[[i]],
              Tx = CMPS_hamby44_results$Tx[[i]],
              include = "nseg",
              outlength = CMPS_hamby44_results$outlength[[i]] # new added
            )
          }
        )
      
      # cat(cb.idx, "- 3; ")
      
      tmp.comp <- tmp.comp %>%
        mutate(
          cmps_score = sapply(tmp.comp$cmps, function(x)
            x$CMPS.score),
          cmps_nseg = sapply(tmp.comp$cmps, function(x)
            x$nseg),
          cmps_score_scaled = cmps_score / cmps_nseg
        )
      
      cmps.table <-
        tmp.comp %>% dplyr::select(land1, land2, landidx1, landidx2, 
                                   cmps_score, cmps_score_scaled, cmps_nseg)
      
      tibble(
        bullet1 = b.cb[, cb.idx][1],
        bullet2 = b.cb[, cb.idx][2],
        cmps.table = list(cmps.table)
      )
    }, mc.cores = detectCores())
  })
  
  # user  system elapsed
  # 0.02    0.00  147.06
  
  
  hamby44.cmps <- do.call(rbind, tmp.44.list)
  
  # stopCluster(cl)
  
  hamby44.cmps$cmps.table.m <-
    lapply(hamby44.cmps$cmps.table, function(tt) {
      tt.idx <- tt %>% rowid_to_column %>%
        filter(land1 %in% idf.scan_id | 
                 land2 %in% idf.scan_id) %>% pull(rowid)
      tt[tt.idx, "cmps_score"] <- NA
      tt[tt.idx, "cmps_score_scaled"] <- NA
      tt
    })
  
  
  # after removing tank-rashed data
  hamby44.cmps <- hamby44.cmps %>% mutate(
    # metric_list = cmps.table.m %>% purrr::map(.f = function(t) {cmps_metrics_helper(t)}) 
    metric_list = cmps.table.m %>% purrr::map(.f = cmps_metrics_helper) 
  )
  
  hamby44.cmps <- bind_cols(hamby44.cmps %>% dplyr::select(-metric_list), 
                             bind_rows(hamby44.cmps$metric_list)) 
  
  # add groundtruth
  hamby44.cmps$type_truth <- tmp.table.44$type_truth
  hamby44.cmps$type <- tmp.table.44$type
  
  # save result in a list
  CMPS_hamby44_results$cmps.table[[i]] <- hamby44.cmps
  
}
saveRDS(CMPS_hamby44_results, "./code/saved_rds/h44-npeaks_aug30.rds")

#############################
# plotting
table.cmps.diff <- lapply(CMPS_hamby44_results$cmps.table, function(tt) {
  type_truth <- tt$type_truth
  # r1 <- compute_var_ratio_anova(tt$cmps.diff, tt$type_truth, MS=FALSE)
  # r2 <- compute_var_ratio_anova(tt$cmps.diff_scaled, tt$type_truth, MS=FALSE)
  # c(r1, r2)
  tt %>% select(cmps.diff:cmps.maxbar_scaled) %>% 
    apply(2, compute_var_ratio_anova, label=type_truth, MS=FALSE)
}) %>% do.call(rbind, .) %>% as.data.frame()

# seg_length
table.cmps.diff$seg_length <- unlist(CMPS_hamby44_results$seg_length)

table.cmps.diff %>%
  pivot_longer(cols = c(cmps.diff:cmps.maxbar_scaled),
               names_to = "score") %>%
  ggplot(aes(x = seg_length, y = value)) +
  geom_line(aes(color = score))

table.cmps.diff$seg_length

metric_plot_helper(CMPS_hamby44_results$cmps.table[[5]], "cmps.diff")
metric_plot_helper(CMPS_hamby44_results$cmps.table[[5]], "cmps.diff_scaled", scaled = TRUE,
                   breaks = seq(0, 1, 0.1), binwidth = 0.05)
metric_plot_helper(CMPS_hamby44_results$cmps.table[[3]], "cmps.diff_scaled", scaled = TRUE)

metric_plot_helper(CMPS_hamby44_results$cmps.table[[7]], "cmps.diff.med_scaled", scaled = TRUE)
metric_plot_helper(CMPS_hamby44_results$cmps.table[[7]], "cmps.diff.med", scaled = FALSE)

metric_plot_helper(CMPS_hamby44_results$cmps.table[[8]], "cmps.max_scaled", scaled = TRUE,
                   breaks = seq(0, 1, 0.1))
metric_plot_helper(CMPS_hamby44_results$cmps.table[[8]], "cmps.max", scaled = FALSE)

metric_plot_helper(CMPS_hamby44_results$cmps.table[[8]], "cmps.diff.med_scaled", scaled = TRUE,
                   breaks = seq(0, 1, 0.1))
metric_plot_helper(CMPS_hamby44_results$cmps.table[[8]], "cmps.diff.med", scaled = FALSE)

CMPS_hamby44_results$cmps.table[[8]] %>% filter(cmps.diff.med == 1, type_truth == "KM") %>% pull(cmps.table) %>% 
  first() %>% 
  ggplot() +
  geom_tile(aes(x = landidx1, y=landidx2, fill = cmps_score)) + 
  scale_fill_gradient2("CMPS score", low = "gray80", high = "darkorange")


########
# npeak
table.cmps.diff$npeak <- unlist(CMPS_hamby44_results$signame)
table.cmps.diff %>% 
  mutate(
    npeak = factor(
      npeak,
      levels = npeak[order(cmps.diff_scaled, decreasing = FALSE)])
  ) %>% 
  pivot_longer(cols = c(cmps.diff:cmps.maxbar_scaled),
               names_to = "score") %>% 
  ggplot(aes(x=npeak, y=value)) +
  geom_bar(aes(fill = score),position = "dodge", stat = "identity") + 
  ylab("variance ratio") +
  coord_flip() + 
  theme_bw()

metric_plot_helper(CMPS_hamby44_results$cmps.table[[7]], 
                   "cmps.diff", scaled = F)



#################################################
# generate plots
for (i in 1:N) {
  hamby44.cmps <- CMPS_hamby44_results$cmps.table[[i]]
  
  hamby44.plot.list <- list()
  titlee <- CMPS_hamby44_results$titlee[[i]]
  
  hamby44.plot.list[[1]] <- hamby44.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.diff,
                       fill = as.factor(type_truth)), binwidth = 1) +
    labs(
      x = expression(bar(CMPS)[diff]),
      fill = "Comparison Type",
      # title = expression(paste(
      #   "Hamby44 - ", bar(CMPS)[max], " Distribution"
      # )),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = seq(0, 24, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  hamby44.plot.list[[2]] <- hamby44.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.diff_scaled,
                       fill = as.factor(type_truth)), binwidth = 0.04) +
    labs(
      x = expression(bar(CMPS^"*")[diff]),
      fill = "Comparison Type",
      # title = expression(paste(
      #   "Hamby44 - ", bar(CMPS)[max], " Distribution"
      # )),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.05)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  plot <- ggarrange(plotlist = hamby44.plot.list,
                    nrow = 1,
                    ncol = 2,
                    common.legend = TRUE, legend = "bottom")
  plot <- annotate_figure(plot, 
                          top = text_grob(com.title44))
  CMPS_hamby44_results$plot[[i]] <- plot
}

# saveRDS(CMPS_hamby44_results, "./code/saved_rds/h44-neapk.rds")


CMPS_hamby44_results <- readRDS("~/Research/CMPSpaper/code/saved_rds/h44-neapk.rds")
CMPS_hamby44_results$plot[[9]]
CMPS_hamby44_results$plot[[3]]
CMPS_hamby44_results$plot[[5]]

# compute variance ratio and generate plots
table.cmps.diff.44 <- lapply(CMPS_hamby44_results$cmps.table, function(tt) {
  r1 <- compute_var_ratio_anova(tt$cmps.diff, tt$type_truth, MS=FALSE)
  r2 <- compute_var_ratio_anova(tt$cmps.diff_scaled, tt$type_truth, MS=FALSE)
  c(r1, r2)
}) %>% do.call(rbind, .) %>% as.data.frame()

colnames(table.cmps.diff.44) <- c("cmps.diff", "cmps.diff_scaled")

# npeak
table.cmps.diff.44$npeak <- unlist(CMPS_hamby44_results$signame)
table.cmps.diff.44 %>% 
  mutate(
    npeak = factor(
      npeak,
      levels = npeak[order(cmps.diff_scaled, decreasing = TRUE)])
  ) %>% 
  pivot_longer(cols = c(cmps.diff, cmps.diff_scaled),
               names_to = "score") %>% 
  ggplot(aes(x=npeak, y=value)) +
  geom_bar(aes(fill = score),position = "dodge", stat = "identity")

# seg_length
table.cmps.diff.44$seg_length <- unlist(CMPS_hamby44_results$seg_length)

table.cmps.diff.44 %>%
  pivot_longer(cols = c(cmps.diff, cmps.diff_scaled),
               names_to = "score") %>%
  ggplot(aes(x = seg_length, y = value)) +
  geom_line(aes(color = score))

order(table.cmps.diff.44$cmps.diff_scaled, decreasing = TRUE)
table.cmps.diff.44$cmps.diff
table.cmps.diff.44$cmps.diff_scaled







