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

data_path <- "./data-csv/hamby44/"


## Data Processing

#### load the data from reproducible/bullet_signatures_etc/
b44 <- read_rds("./bullet_signatures_etc/BulletSignatures44.rds")
b44$sig_rf <- b44$sigs25_531

rf.model <- readRDS("./csafe_rf2.rds")

idf.idx <- c(23,35,39,41,102,126,149)
idf.scan_id <- b44 %>% slice(idf.idx) %>% .$scan_id

# generate bullet_id
bulletid.tb <- b44 %>% select(scan_id)
bulletid.tb <- bulletid.tb %>% mutate(
  barrel_id = sapply(strsplit(scan_id, "-"), "[[", 1),
  bullet = sapply(strsplit(scan_id, "-"), "[[", 2),
  land_id = sapply(strsplit(scan_id, "-"), "[[", 3),
  bullet_id = paste(barrel_id, bullet)
)

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

### set up parameters
span <- 25
N <- 1
CMPS_hamby44_results <- list()
CMPS_hamby44_results$span1 <- as.list(rep(span/100, N))
CMPS_hamby44_results$signame <- list("sig_rf")

CMPS_hamby44_results$titlee <- list("rf_result")
CMPS_hamby44_results$filename <- list("hamby44_rf_result")

CMPS_hamby44_results$rf.table <- list()
CMPS_hamby44_results$plot <- list()

com.title44 <- expression(paste(
  "Hamby 44 - Distribution of various CMPS metrics"
))

# start of the loop
i <- 1

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
#   })

# process of removing outliers
tt <- lapply(b44[[CMPS_hamby44_results$signame[[i]]]], function(x) {
  x$sig
}) %>% unlist()
qtt <- quantile(tt, na.rm = TRUE)
multi.iqr <- 3
outrange <- c(qtt[2] - multi.iqr * (qtt[4] - qtt[2]),
              qtt[4] + multi.iqr * (qtt[4] - qtt[2]))

# outrange <- c(-10000, 10000) # without excluding the outliers

b44 <- b44 %>% mutate(sigs_main = purrr::map(
  .x = .data[[CMPS_hamby44_results$signame[[i]]]],
  .f = function(x) {
    tmp <- x$sig
    tmp[tmp < outrange[1] | tmp > outrange[2]] <- NA
    list(sig = tmp)
  }
))

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
    
    # compute aligned
    tmp.comp <- tmp.comp %>% mutate(
      aligned = purrr::map2(.x = land1, .y = land2, 
                            .f = function(xx, yy) {
                              land1 <- b44$sigs_main[b44$scan_id == xx][[1]]
                              land2 <- b44$sigs_main[b44$scan_id == yy][[1]]
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
                               .f = my_extract_feature_all, resolution = 1.29), 
        legacy_features = purrr::map(striae, 
                                     extract_features_all_legacy, resolution = 1.29))
    # tidyr::unnest(features) # change: instead of legacy_feature
    
    
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

hamby44.rf <- do.call(rbind, tmp.44.list)

# remove outliers/damaged data
hamby44.rf$rf.table.m <-
  lapply(hamby44.rf$rf.table, function(tt) {
    tt.idx <- tt %>% rowid_to_column %>%
      filter(land1 %in% idf.scan_id | land2 %in% idf.scan_id) %>% pull(rowid)
    tt[tt.idx, "rf_score"] <- NA
    
    tt
  })

# after removing tank-rashed data
hamby44.rf <- hamby44.rf %>% mutate(
  metric_list_scaled = rf.table.m %>% purrr::map(.f = function(tt) {
    compute_score_metrics(tt$landidx1, tt$landidx2, tt$rf_score,
                          out_names = c("rf.diff", "rf.diff.med", 
                                        "rf.max", "rf.maxbar"))
  })
) %>% tidyr::unnest(metric_list_scaled)


# add type and type_truth
hamby44.rf$type <- tmp.table.44$type
hamby44.rf$type_truth <- tmp.table.44$type_truth

hamby44.csv <- hamby44.rf %>% select(-rf.table, -rf.table.m)

write.csv(
  hamby44.csv %>% as.data.frame(),
  file = paste(data_path, "hamby44_rf_results", ".csv", sep = "")
)


