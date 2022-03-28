###############################
# Note: please read reproducible-readme.md first before reproducing the results
# Please set the working directory properly so that func_collection.R and 
# the data files mentioned in reproducible-readme.md are available
###############################

# set the working directory
# setwd("~/your-path/supplementary-files")

## Load the required packages
library(tidyverse)
if(!require(bulletxtrctr)) {
  devtools::install_github("heike/bulletxtrctr", ref = "develop")
  library(bulletxtrctr)
}
library(x3ptools)
library(cmpsR)
library(ggpubr)
library(parallel)

# source("func_collection.R")

data_path <- "./data-csv/hamby44/"

## Data Processing

#### load the data from reproducible/bullet_signatures_etc/
b44 <- read_rds("./bullet_signatures_etc/BulletSignatures44.rds")

#### obtain all comparisons
bulletid.tb <- b44 %>% select(scan_id)
bulletid.tb <- bulletid.tb %>% mutate(
  barrel_id = sapply(strsplit(scan_id, "-"), "[[", 1),
  bullet = sapply(strsplit(scan_id, "-"), "[[", 2),
  land_id = sapply(strsplit(scan_id, "-"), "[[", 3),
  bullet_id = paste(barrel_id, bullet)
)
b44$bullet_id <- bulletid.tb %>% pull(bullet_id)
b.cb <- unique(b44$bullet_id) %>% utils::combn(m = 2)
p <- 1:595

#### import indices of damaged scans
###### These indices are identified by using a Shiny App called bulletinspectR
idf.idx <- c(23,35,39,41,102,126,149)
idf.scan_id <- b44 %>% slice(idf.idx) %>% .$scan_id

#### include the ground truth for Hamby44
key.table <- tibble(
  bullet.id = c("K",
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
                "S"),
  barrel.id = c(1, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 9, 10)
)

tmp.table.44 <- tibble(bullet1 = b.cb[1,],
                       bullet2 = b.cb[2,])

tmp.table.44$b1t <-
  tmp.table.44$bullet1 %>% str_extract("'[A-Z]'") %>% str_extract("[A-Z]")
tmp.table.44$b2t <-
  tmp.table.44$bullet2 %>% str_extract("'[A-Z]'") %>% str_extract("[A-Z]")

tmp.table.44$b1b <- tmp.table.44 %>% 
  left_join(key.table, by = c("b1t" = "bullet.id")) %>% pull(barrel.id)
tmp.table.44$b2b <- tmp.table.44 %>% 
  left_join(key.table, by = c("b2t" = "bullet.id")) %>% pull(barrel.id)

####### KM - come from the same barrel; KNM - otherwise
tmp.table.44 <- tmp.table.44 %>% mutate(
  b1b = if_else(
    is.na(b1b),
    bullet1 %>% str_extract("Barrel \\d+") %>% 
      str_extract("\\d+") %>% as.numeric(),
    b1b
  ),
  b2b = if_else(
    is.na(b2b),
    bullet2 %>% str_extract("Barrel \\d+") %>% 
      str_extract("\\d+") %>% as.numeric(),
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

## set up parameters
N <- 8 # the number of different sets of parameters
CMPS_hamby44_results <- list() # a container for everything

#### setup span1
CMPS_hamby44_results$span1 <- as.list(rep(0.25, N))

#### setup signature name
CMPS_hamby44_results$signame <-
  as.list(rep("sigs25_531", N))

#### setup neapks.set
CMPS_hamby44_results$npeaks_set <-
  lapply(1:N, function(t) c(5,3,1))

#### setup seg_length
CMPS_hamby44_results$seg_length <- as.list(c(30, 45, 61, 90, 122, 150, 180, 210))

#### setup outlength
CMPS_hamby44_results$outlength <- vector(mode = "list", length = N)

#### setup Tx
CMPS_hamby44_results$Tx <- as.list(rep(30, N))

#### create a plot title and a file name for each set of parameters
CMPS_hamby44_results$titlee <- list()
CMPS_hamby44_results$filename <- list()
for (i in 1:N) {
  CMPS_hamby44_results$titlee[[i]] <-
    paste0(
      "npeaks_set=c(",
      paste(CMPS_hamby44_results$npeaks_set[[i]], collapse = ","),
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
      "hamby44_segment",
      CMPS_hamby44_results$span1[[i]]*100,
      paste(CMPS_hamby44_results$npeaks_set[[i]], collapse = "-"),
      CMPS_hamby44_results$seg_length[[i]],
      CMPS_hamby44_results$Tx[[i]],
      sep = "_"
    )
}

#### container for the cmps results and corresponding plots
CMPS_hamby44_results$cmps.table <- list()
CMPS_hamby44_results$plot <- list()

#### a global title for Hamby44 plots
com.title44 <- expression(paste(
  "Hamby 44 - ", CMPS[max], " and ", bar(CMPS)[max], " Distribution"
))

## Compute CMPS scores

#### for each set of parameters, compute the CMPS scores
for (i in 1:N) {
  
  ###### code for generating bullet signatures; commented out
  # b44[[CMPS_hamby44_results$signame[[i]]]] <- purrr::map2(
  #   .x = b44.full$ccdata,
  #   .y = b44.full$grooves,
  #   .f = function(x, y) {
  #     cc_get_signature(
  #       ccdata = x,
  #       grooves = y,
  #       span1 = CMPS_hamby44_results$span1[[i]],
  #       span2 = 0.03
  #     )
  #   }
  # )
  
  #### remove outliers in bullet signatures
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
  
  cl <- parallel::makeCluster(detectCores())
  par.setup <- parLapply(cl, 1:length(cl),
                         function(xx) {
                           library(tidyverse)
                           library(cmpsR)
                         })
  clusterExport(cl, list('b44', 'b.cb', 'CMPS_hamby44_results','i'))
  
  #### compute CMPS scores with parallel computing
  system.time({
    tmp.44.list <- parLapply(cl, p, function(cb.idx) {
    # tmp.44.list <- mclapply(p, function(cb.idx) {
      
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
              npeaks_set = CMPS_hamby44_results$npeaks_set[[i]],
              seg_length = CMPS_hamby44_results$seg_length[[i]],
              Tx = CMPS_hamby44_results$Tx[[i]],
              include = "nseg",
              outlength = CMPS_hamby44_results$outlength[[i]] 
            )
          }
        )
      
      # cat(cb.idx, "- 3; ")
      
      tmp.comp <- tmp.comp %>%
        mutate(
          cmps_score = sapply(tmp.comp$cmps, function(x)
            x$CMPS_score),
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
    # }, mc.cores = detectCores())
    })
  })
  
  # user  system elapsed
  # 0.02    0.00  147.06
  hamby44.cmps <- do.call(rbind, tmp.44.list)
  
  stopCluster(cl)
  
  ## compute CMPS_{max} and \bar{CMPS_{max}} for the bullet comparison
  # set NAs for tank-rashed data
  hamby44.cmps$cmps.table.m <-
    lapply(hamby44.cmps$cmps.table, function(tt) {
      tt.idx <- tt %>% rowid_to_column %>%
        filter(land1 %in% idf.scan_id |
                 land2 %in% idf.scan_id) %>% pull(rowid)
      tt[tt.idx, "cmps_score"] <- NA
      tt[tt.idx, "cmps_score_scaled"] <- NA
      
      tt
    })
  
  # after removing tank-rashed data, recompute CMPS_{max} and \bar{CMPS_{max}}
  hamby44.cmps <- hamby44.cmps %>% mutate(
    metric_list = cmps.table.m %>% purrr::map(.f = function(tt) {
      compute_score_metrics(tt$landidx1, tt$landidx2, tt$cmps_score,
                            out_names = c("cmps.diff", "cmps.diff.med", "cmps.max", "cmps.maxbar"))
    }),
    metric_list_scaled = cmps.table.m %>% purrr::map(.f = function(tt) {
      compute_score_metrics(tt$landidx1, tt$landidx2, tt$cmps_score_scaled,
                            out_names = c("cmps.diff_scaled", "cmps.diff.med_scaled", 
                                          "cmps.max_scaled", "cmps.maxbar_scaled"))
    })
  )
  
  hamby44.cmps <- hamby44.cmps %>% tidyr::unnest(c(metric_list, metric_list_scaled))
  
  # add ground truth and type
  hamby44.cmps$type_truth <- tmp.table.44$type_truth
  hamby44.cmps$type <- tmp.table.44$type
  
  # save result in a list
  CMPS_hamby44_results$cmps.table[[i]] <- hamby44.cmps
  
}

#### code used to save results as csv files
for(i in 1:N){
  write.csv(
    CMPS_hamby44_results$cmps.table[[i]] %>% select(-cmps.table, -cmps.table.m) %>% as.data.frame(),
    file = paste(data_path, CMPS_hamby44_results$filename[[i]], ".csv", sep = ""))
}

## Generate plots
# for (i in 1:N) {
#   hamby44.cmps <- CMPS_hamby44_results$cmps.table[[i]]
#   
#   hamby44.plot.list <- list()
#   titlee <- CMPS_hamby44_results$titlee[[i]]
#   
#   # generate plot for CMPS_{max}
#   hamby44.plot.list[[1]] <- hamby44.cmps %>% ggplot() +
#     geom_histogram(aes(x = cmps.max.m,
#                        fill = as.factor(type_truth)), binwidth = 1) +
#     labs(
#       fill = "Comparison Type",
#       x = expression(CMPS[max]),
#       subtitle = titlee
#     ) +
#     scale_x_continuous(breaks = seq(0, 27, 1)) +
#     theme_bw() +
#     theme(panel.grid.minor = element_blank()) +
#     font("x.text", size = 6)
#   
#   # generate plot for \bar{CMPS_{max}}
#   hamby44.plot.list[[2]] <- hamby44.cmps %>% ggplot() +
#     geom_histogram(aes(x = cmps.maxbar.m,
#                        fill = as.factor(type_truth)), binwidth = 1) +
#     labs(
#       x = expression(bar(CMPS)[max]),
#       fill = "Comparison Type",
#       subtitle = titlee
#     ) +
#     scale_x_continuous(breaks = seq(0, 24, 1)) +
#     theme_bw() +
#     theme(panel.grid.minor = element_blank())
#   
#   # combine the two plots
#   plot <- ggarrange(plotlist = hamby44.plot.list,
#                     nrow = 1,
#                     ncol = 2,
#                     common.legend = TRUE, legend = "bottom")
#   plot <- annotate_figure(plot, 
#                           top = text_grob(com.title44))
#   CMPS_hamby44_results$plot[[i]] <- plot
# }
