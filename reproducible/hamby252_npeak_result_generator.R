##########################
# Note: This R script uses rds files in reproducible/bullet_signatures_etc/
#       These rds files contain processed bullet signatures used in the paper
#       The original data contains x3p objects of large size, so they are not
#       uploaded to the Github repo
##########################

## Load the required packages
library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(CMPS)
library(ggpubr)
library(parallel)

source("code/func_collection.R")

#### code used to load the original data; commented out
# b252.full <- 
#   readRDS("~/Research/CMPSpaper/preconsideration/bullets_hamby252_Sep9.rds")

## Data Processing

#### load the data from reproducible/bullet_signatures_etc/
b252 <- read_rds("./reproducible/bullet_signatures_etc/BulletSignatures252.rds")

#### obtain all comparisons
b.cb <- unique(b252$bullet_id) %>% utils::combn(m = 2)
p <- 1:595

#### obtained ID for the damaged bullets;
###### These IDs are identified by using a Shiny App called bulletinspectR
tr.id <- c("Br3-1-4", "Ukn-B-2", "Ukn-Q-4", "Br6-2-1", "Br9-2-4")

#### Include the ground truth for Hamby252
tmp.tibble <- tibble(bullet1 = b.cb[1,], 
                     bullet2 = b.cb[2,])

tmp.tibble <- tmp.tibble %>% mutate(
  type = if_else(
    bullet1 %>% str_extract("Br\\d+") ==
      bullet2 %>% str_extract("Br\\d+"),
    "KM",
    "KNM"
  ),
  type = if_else(is.na(type), "Ukn", type)
)

###### StudyInfo.xlsx provides the ground truth for Hamby252 
studyinfo <-
  readxl::read_xlsx(
    "./reproducible/bullet_signatures_etc/StudyInfo.xlsx",
    sheet = 3,
    skip = 1
  )

keyss <- studyinfo$`Specimen ID`[-length(studyinfo$`Specimen ID`)]
keyss <- tibble(specimen = keyss)
keyss <- keyss %>% mutate(
  barrel = specimen %>% str_split(" "),
  barrel.id = purrr::map_chr(barrel, function(bb) {
    bb[1]
  }),
  bullet.id = purrr::map_chr(barrel, function(bb) {
    bb[3]
  })
) %>% select(-barrel)

###### obtain the matching keys
keyss <- keyss %>% filter(!(bullet.id %>% str_detect("\\d+")))

tmp.ukn <-
  tibble(
    col1 = tmp.tibble$bullet1 %>% str_extract(" [A-Z]") %>% str_trim(),
    col2 = tmp.tibble$bullet2 %>% str_extract(" [A-Z]") %>% str_trim()
  )

###### left join the ground truth
tmp.tibble$b1t <-
  tmp.ukn %>% left_join(keyss, by = c("col1" = "bullet.id")) %>%
  pull(barrel.id) %>% str_remove("Brl")
tmp.tibble$b2t <-
  tmp.ukn %>% left_join(keyss, by = c("col2" = "bullet.id")) %>%
  pull(barrel.id) %>% str_remove("Brl")

####### KM - come from the same barrel; KNM - otherwise
tmp.tibble <- tmp.tibble %>% mutate(
  b1t = if_else(
    is.na(b1t),
    bullet1 %>% str_extract("Br\\d+") %>% str_remove("Br"),
    b1t
  ),
  b2t = if_else(
    is.na(b2t),
    bullet2 %>% str_extract("Br\\d+") %>% str_remove("Br"),
    b2t
  ),
  b1t = as.numeric(b1t),
  b2t = as.numeric(b2t),
  type_truth = if_else(b1t == b2t, "KM", "KNM")
)

## set up parameters
N <- 8 # the number of different sets of parameters
CMPS_hamby252_results <- list() # a container for everything

#### setup span1
CMPS_hamby252_results$span1 <- as.list(rep(0.25, N))

#### setup signature name
CMPS_hamby252_results$signame <- as.list(rep("sigs25", N))
# list("sigs75", "sigs25", "sigs25_1062", "sigs15")

#### setup npeaks.set
CMPS_hamby252_results$npeaks.set <- 
  list(c(5),
       c(5,2),
       c(5,3,1),
       c(10,6,2),
       c(6,4,2,1),
       c(10,7,4,2),
       c(10,6,4,2,1),
       c(10,8,6,4,2,1))

#### setup seg_length
CMPS_hamby252_results$seg_length <- 
  # as.list(c(25, 50, 75, 100, 125, 150, 175, 200))
  as.list(rep(50, N))

#### setup Tx
CMPS_hamby252_results$Tx <- as.list(rep(25, N))

#### setup outlength
CMPS_hamby252_results$outlength <- vector(mode = "list", length = N)
CMPS_hamby252_results$outlength[[7]] <- c(60, 90, 120, 150, 180)
CMPS_hamby252_results$outlength[[8]] <- c(50,75,100,125,150,200)

#### create a plot title and a file name for each set of parameters
CMPS_hamby252_results$titlee <- list()
CMPS_hamby252_results$filename <- list()
for (i in 1:N) {
  CMPS_hamby252_results$titlee[[i]] <-
    paste0(
      "npeaks.set=c(",
      paste(CMPS_hamby252_results$npeaks.set[[i]], collapse = ","),
      ")",
      ", len=",
      CMPS_hamby252_results$seg_length[[i]],
      ", Tx=",
      CMPS_hamby252_results$Tx[[i]],
      ", \nspan1=",
      CMPS_hamby252_results$span1[[i]],
      ", span2=0.03"
    )
  CMPS_hamby252_results$filename[[i]] <- 
    paste(
      "hamby252_npeak",
      CMPS_hamby252_results$span1[[i]]*100,
      paste(CMPS_hamby252_results$npeaks.set[[i]], collapse = "-"),
      CMPS_hamby252_results$seg_length[[i]],
      CMPS_hamby252_results$Tx[[i]],
      sep = "_"
    )
}

#### container for the cmps results and corresponding plots
CMPS_hamby252_results$cmps.table <- list()
CMPS_hamby252_results$plot <- list()

#### a global title for Hamby252 plots
com.title252 <- expression(paste(
  "Hamby 252 - ", CMPS[max], " and ", bar(CMPS)[max], " Distribution"
))

## Compute CMPS scores

#### for each set of parameters, compute the CMPS scores
for (i in 1:N) {
  
  ###### code for generating bullet signatures; commented out
  # b252[[CMPS_hamby252_results$signame[[i]]]] <- purrr::map2(
  #   .x = b252.full$ccdata,
  #   .y = b252.full$grooves,
  #   .f = function(x, y) {
  #     cc_get_signature(
  #       ccdata = x,
  #       grooves = y,
  #       span1 = CMPS_hamby252_results$span1[[i]],
  #       span2 = 0.03
  #     )
  #   }
  # )
  
  #### remove outliers in bullet signatures
  tt <-
    lapply(b252[[CMPS_hamby252_results$signame[[i]]]], function(x) {
      x$sig
    }) %>% unlist()
  qtt <- quantile(tt, na.rm = TRUE)
  multi.iqr <- 3
  outrange <-
    c(
      qtt[2] - multi.iqr * (qtt[4] - qtt[2]), 
      qtt[4] + multi.iqr * (qtt[4] - qtt[2])
    )
  b252 <- b252 %>% mutate(sigs_main = purrr::map(
    .x = .data[[CMPS_hamby252_results$signame[[i]]]],
    .f = function(x) {
      tmp <- x$sig
      tmp[tmp < outrange[1] | tmp > outrange[2]] <- NA
      list(sig = tmp)
    }
  ))
  
  
  #### compute CMPS scores with parallel computing
  system.time({
    tmp.252.list <- mclapply(p, function(cb.idx) { # for each pairwise comparison
      
      # identify two bullets with their IDs
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
      
      # compute CMPS scores for each signature comparison
      tmp.comp$cmps <-
        purrr::map2(
          .x = tmp.comp$land1,
          .y = tmp.comp$land2,
          .f = function(xx, yy) {
            sig1 <- b252$sigs_main[b252$scan_id == xx][[1]]$sig
            sig2 <- b252$sigs_main[b252$scan_id == yy][[1]]$sig
            
            extract_feature_cmps(
              sig1,
              sig2,
              npeaks.set = CMPS_hamby252_results$npeaks.set[[i]],
              seg_length = CMPS_hamby252_results$seg_length[[i]],
              Tx = CMPS_hamby252_results$Tx[[i]],
              include = "nseg",
              outlength = CMPS_hamby252_results$outlength[[i]]
            )
          }
        )
      
      # extract CMPS scores
      tmp.comp <- tmp.comp %>%
        mutate(
          cmps_score = sapply(tmp.comp$cmps, function(x)
            x$CMPS.score),
          cmps_nseg = sapply(tmp.comp$cmps, function(x)
            x$nseg),
          cmps_score_scaled = cmps_score / cmps_nseg
        )
      
      # summarize results and return
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
  # 0.02    0.00  126.06
  
  hamby252.cmps <- do.call(rbind, tmp.252.list)
  
  # compute CMPS_{max} and \bar{CMPS_{max}} for the bullet comparison
  # set NAs for tank-rashed data
  hamby252.cmps$cmps.table.m <-
    lapply(hamby252.cmps$cmps.table, function(tt) {
      tt.idx <- tt %>% rowid_to_column %>%
        filter(land1 %in% tr.id | land2 %in% tr.id) %>% pull(rowid)
      tt[tt.idx, "cmps_score"] <- NA
      tt[tt.idx, "cmps_score_scaled"] <- NA
      
      tt
    })
  
  # after removing tank-rashed data, recompute CMPS_{max} and \bar{CMPS_{max}}
  hamby252.cmps <- hamby252.cmps %>% mutate(
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
  
  hamby252.cmps <- hamby252.cmps %>% tidyr::unnest(c(metric_list, metric_list_scaled))
  
  # add ground truth and type
  hamby252.cmps$type <- tmp.tibble$type
  hamby252.cmps$type_truth <- tmp.tibble$type_truth
  
  # save result in a list
  CMPS_hamby252_results$cmps.table[[i]] <- hamby252.cmps
  
}

#### code used to save results as csv files used in the paper; commented out
data_path <- "~/Research/CMPSpaper/CMPSpaper_writing/data/hamby252/"
for(i in 1:N){
  write.csv(
    CMPS_hamby252_results$cmps.table[[i]] %>%
      select(-cmps.table, -cmps.table.m) %>% as.data.frame(),
    file = paste(data_path, CMPS_hamby252_results$filename[[i]], ".csv", sep = ""))
}

## Generate plots
for (i in 1:N) {
  hamby252.cmps <- CMPS_hamby252_results$cmps.table[[i]]
  
  hamby252.plot.list <- list()
  titlee <- CMPS_hamby252_results$titlee[[i]]
  
  # generate plot for CMPS_{max}
  hamby252.plot.list[[1]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.max.m,
                       fill = as.factor(type_truth)), binwidth = 1) +
    labs(
      fill = "Comparison Type",
      x = expression(CMPS[max]),
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = seq(0, 27, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    font("x.text", size = 6)
  
  # generate plot for \bar{CMPS_{max}}
  hamby252.plot.list[[2]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.maxbar.m,
                       fill = as.factor(type_truth)), binwidth = 1) +
    labs(
      x = expression(bar(CMPS)[max]),
      fill = "Comparison Type",
      subtitle = titlee
    ) +
    scale_x_continuous(breaks = seq(0, 24, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  # combine the two plots
  plot <- ggarrange(plotlist = hamby252.plot.list,
                    nrow = 1,
                    ncol = 2,
                    common.legend = TRUE, legend = "bottom")
  plot <- annotate_figure(plot, 
                          top = text_grob(com.title252))
  CMPS_hamby252_results$plot[[i]] <- plot
}
