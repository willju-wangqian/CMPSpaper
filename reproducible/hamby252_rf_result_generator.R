###############################
# Note: read reproducible-readme.md first before reproducing the results
# Please set the working directory properly so that func_collection.R and 
# the data files mentioned in reproducible-readme.md are available
###############################

# set the working directory
# setwd("~/your-path/supplementary-files")

## Load the required packages
library(tidyverse)
if(!require(bulletxtrctr)) {
  devtools::install_github("heike/bulletxtrctr")
  library(bulletxtrctr)
}
library(x3ptools)
library(cmpsR)
library(ggpubr)
library(parallel)
library(assertthat)
require(randomForest)

source("func_collection.R")

data_path <- "./data-csv/hamby252/"

#### load the data from reproducible/bullet_signatures_etc/
b252 <- read_rds("./bullet_signatures_etc/BulletSignatures252.rds")
b252$sig_rf <- b252$sigs25

rf.model <- readRDS("./csafe_rf2.rds")

# generate bullet_id
bulletid.tb <- b252 %>% select(scan_id)
bulletid.tb <- bulletid.tb %>% mutate(
  barrel_id = sapply(strsplit(scan_id, "-"), "[[", 1),
  bullet = sapply(strsplit(scan_id, "-"), "[[", 2),
  land_id = sapply(strsplit(scan_id, "-"), "[[", 3),
  bullet_id = paste(barrel_id, bullet)
)

# get bullet id
b252$bullet_id <- bulletid.tb %>% pull(bullet_id)

# get all comparisons
b.cb <- unique(b252$bullet_id) %>% utils::combn(m = 2)
p <- 1:595

# remove damaged data
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

tr.id <- c("Br3-1-4", "Ukn-B-2", "Ukn-Q-4", "Br6-2-1", "Br9-2-4")

######################################
# add ground-truth for hamby 252
studyinfo <-
  readxl::read_xlsx(
    "./bullet_signatures_etc/StudyInfo.xlsx",
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

# obtain the matching keys
keyss <- keyss %>% filter(!(bullet.id %>% str_detect("\\d+")))

tmp.ukn <-
  tibble(
    col1 = tmp.tibble$bullet1 %>% str_extract(" [A-Z]") %>% str_trim(),
    col2 = tmp.tibble$bullet2 %>% str_extract(" [A-Z]") %>% str_trim()
  )

# left join the ground truth
tmp.tibble$b1t <-
  tmp.ukn %>% left_join(keyss, by = c("col1" = "bullet.id")) %>%
  pull(barrel.id) %>% str_remove("Brl")
tmp.tibble$b2t <-
  tmp.ukn %>% left_join(keyss, by = c("col2" = "bullet.id")) %>%
  pull(barrel.id) %>% str_remove("Brl")

# find KM KNM
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

### set up parameters
span <- 25
N <- 1
CMPS_hamby252_results <- list()
CMPS_hamby252_results$span1 <- as.list(rep(span/100, N))
CMPS_hamby252_results$signame <- list("sig_rf")

CMPS_hamby252_results$titlee <- list("rf_result")
CMPS_hamby252_results$filename <- list("hamby252_rf_result")

CMPS_hamby252_results$rf.table <- list()
CMPS_hamby252_results$plot <- list()

com.title252 <- expression(paste(
  "Hamby 252 - Distribution of various CMPS metrics"
))

# start of the loop
i <- 1

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
#   })


# process of removing outliers
tt <- lapply(b252[[CMPS_hamby252_results$signame[[i]]]], function(x) {
  x$sig
}) %>% unlist()
qtt <- quantile(tt, na.rm = TRUE)
multi.iqr <- 3
outrange <- c(qtt[2] - multi.iqr * (qtt[4] - qtt[2]),
              qtt[4] + multi.iqr * (qtt[4] - qtt[2]))

# outrange <- c(-10000, 10000) # without excluding the outliers

b252 <- b252 %>% mutate(sigs_main = purrr::map(
  .x = .data[[CMPS_hamby252_results$signame[[i]]]],
  .f = function(x) {
    tmp <- x$sig
    tmp[tmp < outrange[1] | tmp > outrange[2]] <- NA
    list(sig = tmp)
  }
))

# start the parallel computing
system.time({
  # tmp.252.list <- parLapply(cl, p, function(cb.idx) {
  tmp.252.list <- mclapply(p, function(cb.idx) {
    
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
    
    # compute aligned
    tmp.comp <- tmp.comp %>% mutate(
      aligned = purrr::map2(.x = land1, .y = land2, 
                            .f = function(xx, yy) {
                              land1 <- b252$sigs_main[b252$scan_id == xx][[1]]
                              land2 <- b252$sigs_main[b252$scan_id == yy][[1]]
                              land1$bullet <- "first-land"
                              land2$bullet <- "second-land"
                              
                              sig_align(land1$sig, land2$sig)
                              # bulletxtrctr::sig_align(land1$sig, land2$sig)
                            }))
    
    # compute the rf_score
    tmp.comp <- tmp.comp %>% 
      mutate(
        striae = aligned %>% purrr::map(.f = sig_cms_max),
        features = purrr::map2(.x = aligned, .y = striae, 
                               .f = my_extract_feature_all, resolution = 1.5625), 
        legacy_features = purrr::map(striae, 
                                     extract_features_all_legacy, resolution = 1.5625))
    # tidyr::unnest(features) # change: instead of legacy_feature
    
    # tmp.comp$rf_score <- predict(bulletxtrctr::rtrees, newdata = tmp.comp, 
    #                             type = "prob")[, 2]
    tmp.comp$rf_score <- predict(rf.model, 
                                 newdata = tmp.comp %>% tidyr::unnest(features) %>% mutate(
                                   abs_lag_mm = abs(lag_mm)
                                 ),
                                 type = "prob")[, 2]
    
    rf.table <-
      tmp.comp %>% dplyr::select(land1, land2, landidx1, landidx2, 
                                 rf_score, features, legacy_features)
    
    tibble(
      bullet1 = b.cb[, cb.idx][1],
      bullet2 = b.cb[, cb.idx][2],
      rf.table = list(rf.table)
    )
    # })
  }, mc.cores = detectCores())
})
# user  system elapsed
# 0.02    0.00  126.06

hamby252.rf <- do.call(rbind, tmp.252.list)

# remove outliers/damaged data
hamby252.rf$rf.table.m <-
  lapply(hamby252.rf$rf.table, function(tt) {
    tt.idx <- tt %>% rowid_to_column %>%
      filter(land1 %in% tr.id | land2 %in% tr.id) %>% pull(rowid)
    tt[tt.idx, "rf_score"] <- NA
    
    tt
  })

# after removing tank-rashed data
hamby252.rf <- hamby252.rf %>% mutate(
  metric_list_scaled = rf.table.m %>% purrr::map(.f = function(tt) {
    compute_score_metrics(tt$landidx1, tt$landidx2, tt$rf_score,
                          out_names = c("rf.diff", "rf.diff.med", 
                                        "rf.max", "rf.maxbar"))
  })
) %>% tidyr::unnest(metric_list_scaled)

# add type and type_truth
hamby252.rf$type <- tmp.tibble$type
hamby252.rf$type_truth <- tmp.tibble$type_truth

hamby252.csv <- hamby252.rf %>% select(-rf.table, -rf.table.m)

write.csv(
  hamby252.csv %>% as.data.frame(),
  file = paste(data_path, "hamby252_rf_results", ".csv", sep = "")
)

