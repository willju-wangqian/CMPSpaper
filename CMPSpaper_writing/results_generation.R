library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(CMPS)
library(ggpubr)
library(parallel)

###########
# from the x3p data, re-run CMPS results for hamby 252
# hamby252_results <-
#   readRDS("C:/Research/Bullet/Bullet_Rcode/hamby252_results.rds")
# hamby252_results <- 
#   readRDS("~/Research/CMPSpaper/preconsideration/hamby252_results.rds")

# b252.full <-
#   readRDS("C:/Research/Bullet/Bullet_Rcode/hamby252/bullets_hamby252_Sep9.rds")
b252.full <- 
  readRDS("~/Research/CMPSpaper/preconsideration/bullets_hamby252_Sep9.rds")

# generate bullet_id
bulletid.tb <- b252.full %>% select(scan_id)
bulletid.tb <- bulletid.tb %>% mutate(
  barrel_id = sapply(strsplit(scan_id, "-"), "[[", 1),
  bullet = sapply(strsplit(scan_id, "-"), "[[", 2),
  land_id = sapply(strsplit(scan_id, "-"), "[[", 3),
  bullet_id = paste(barrel_id, bullet)
)

# remove x3p files
b252 <- b252.full %>% select(-x3p,-crosscut,-grooves,-ccdata)

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

# include tank-rashed data
# [TODO] MAKE THIS A CSV FILE
# t1.df <-
#   # read_rds("C:/Research/Bullet/Bullet_Rcode/hamby252/t1.df.rds")
#   read_rds("~/Research/CMPSpaper/preconsideration/t1.df.rds")
# tr.id <-
#   c(
#     t1.df %>% filter(type == "TR") %>% pull(scan_id) %>% as.character(),
#     "Br6-2-1",
#     "Br9-2-4"
#   )

tr.id <- c("Br3-1-4", "Ukn-B-2", "Ukn-Q-4", "Br6-2-1", "Br9-2-4")

######################################
# add ground-truth for hamby 252
studyinfo <-
  readxl::read_xlsx(
    "~/Research/CMPSpaper/preconsideration/StudyInfo.xlsx",
    sheet = 3,
    skip = 1
  )
# readxl::read_xlsx(
#   "C:/Research/Bullet/Bullet_Rcode/grooves/NRBTDSearchResults/Hamby (2009) Barrel/StudyInfo.xlsx",
#   sheet = 3,
#   skip = 1
# )

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
N <- 4
CMPS_hamby252_results <- list()
CMPS_hamby252_results$span1 <- list(0.75, 0.25, 0.25, 0.15)
CMPS_hamby252_results$signame <-
  list("sigs75", "sigs25", "sigs25_1062", "sigs15")
CMPS_hamby252_results$npeaks.set <-
  list(c(5, 3, 1),
       c(5, 3, 1),
       c(10, 6, 2),
       c(5,3,1))
CMPS_hamby252_results$seg_length <- as.list(rep(50, 4))
CMPS_hamby252_results$Tx <- as.list(rep(25, 4))
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
      "hamby252",
      CMPS_hamby252_results$span1[[i]]*100,
      paste(CMPS_hamby252_results$npeaks.set[[i]], collapse = "-"),
      CMPS_hamby252_results$seg_length[[i]],
      CMPS_hamby252_results$Tx[[i]],
      sep = "_"
    )
}
CMPS_hamby252_results$cmps.table <- list()
CMPS_hamby252_results$plot <- list()

com.title252 <- expression(paste(
  "Hamby 252 - ", CMPS[max], " and ", bar(CMPS)[max], " Distribution"
))


for (i in 1:N) {
  # i <- 1
  
  b252[[CMPS_hamby252_results$signame[[i]]]] <- purrr::map2(
    .x = b252.full$ccdata,
    .y = b252.full$grooves,
    .f = function(x, y) {
      cc_get_signature(
        ccdata = x,
        grooves = y,
        span1 = CMPS_hamby252_results$span1[[i]],
        span2 = 0.03
      )
    }
  )
  
  # process of removing outliers
  tt <-
    lapply(b252[[CMPS_hamby252_results$signame[[i]]]], function(x) {
      x$sig
    }) %>% unlist()
  qtt <- quantile(tt, na.rm = TRUE)
  multi.iqr <- 3
  outrange <-
    c(qtt[2] - multi.iqr * (qtt[4] - qtt[2]), qtt[4] + multi.iqr * (qtt[4] -
                                                                      qtt[2]))
  
  b252 <- b252 %>% mutate(sigs_main = purrr::map(
    .x = .data[[CMPS_hamby252_results$signame[[i]]]],
    .f = function(x) {
      tmp <- x$sig
      tmp[tmp < outrange[1] | tmp > outrange[2]] <- NA
      list(sig = tmp)
    }
  ))
  

  
  ## parallel setup
  # cl <- makeCluster(6)
  # 
  # par.setup <- parLapply(cl, 1:length(cl),
  #                        function(xx) {
  #                          library(tidyverse)
  #                          library(bulletxtrctr)
  #                          library(x3ptools)
  #                          library(CMPS)
  #                        })
  # 
  # clusterExport(cl,
  #               c("b252", "b.cb", "p",
  #                 "CMPS_hamby252_results", "i"),
  #               envir = globalenv())
  
  # start the parallel computing
  system.time({
    # tmp.252.list <- parLapply(cl, p, function(cb.idx) {
    tmp.252.list <- mclapply(p, function(cb.idx) {
      # cat(" ############# \n", "start of ", cb.idx, "\n")
      
      tmp.lands <-
        c(
          b252 %>% filter(bullet_id %in% b.cb[, cb.idx][1]) %>% pull(scan_id),
          b252 %>% filter(bullet_id %in% b.cb[, cb.idx][2]) %>% pull(scan_id)
        )
      tmp.comp <-
        data.frame(expand.grid(land1 = tmp.lands[1:6], land2 = tmp.lands[7:12]),
                   stringsAsFactors = FALSE)
      
      # cat(cb.idx, "- 1; ")
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
              include = "nseg"
            )
          }
        )
      
      # cat(cb.idx, "- 3; ")
      
      tmp.comp <- tmp.comp %>%
        mutate(
          cmps_score = sapply(tmp.comp$cmps, function(x)
            x$CMPS.score),
          cmps_nseg = sapply(tmp.comp$cmps, function(x)
            x$nseg)
        )
      
      # cat(cb.idx, "- 4;\n")
      
      cmps.table <-
        tmp.comp %>% select(land1, land2, cmps_score, cmps_nseg)
      # cp2 <- tmp.comp %>% select(land1, land2, cmps_score, cmps_nseg)
      
      # cat(" end of ", cb.idx, "\n", "############# \n")
      
      tibble(
        bullet1 = b.cb[, cb.idx][1],
        bullet2 = b.cb[, cb.idx][2],
        cmps.table = list(cmps.table)
      )
    # })
    }, mc.cores = detectCores())
  })
  
  # user  system elapsed
  # 0.02    0.00  126.06
  
  
  hamby252.cmps <- do.call(rbind, tmp.252.list)
  
  # stopCluster(cl)
  
  ##################################################
  # before removing tank-rashed data
  hamby252.cmps <- hamby252.cmps %>% mutate(
    cmps.max = cmps.table %>% purrr::map_dbl(
      .f = function(t) {
        max(t$cmps_score)
      }
    ),
    cmps.maxbar = cmps.table %>% purrr::map_dbl(
      .f = function(t) {
        landidx1 <-
          t$land1 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        landidx2 <-
          t$land2 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        compute_average_scores(landidx1, landidx2, t$cmps_score) %>% max()
      }
    )
  )
  
  
  
  hamby252.cmps$cmps.table.m <-
    lapply(hamby252.cmps$cmps.table, function(tt) {
      tt.idx <- tt %>% rowid_to_column %>%
        filter(land1 %in% tr.id | land2 %in% tr.id) %>% pull(rowid)
      tt[tt.idx, "cmps_score"] <- NA
      tt
    })
  
  # after removing tank-rashed data
  hamby252.cmps <- hamby252.cmps %>% mutate(
    cmps.max.m = cmps.table.m %>% purrr::map_dbl(
      .f = function(t) {
        max(t$cmps_score, na.rm = TRUE)
      }
    ),
    cmps.maxbar.m = cmps.table.m %>% purrr::map_dbl(
      .f = function(t) {
        landidx1 <-
          t$land1 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        landidx2 <-
          t$land2 %>% as.character() %>% str_sub(-1, -1) %>% as.numeric()
        compute_average_scores(landidx1, landidx2, t$cmps_score, addNA = TRUE) %>% max()
      }
    )
  )
  
  # add type and type_truth
  hamby252.cmps$type <- tmp.tibble$type
  hamby252.cmps$type_truth <- tmp.tibble$type_truth
  
  ###################################
  # save result in a list
  CMPS_hamby252_results$cmps.table[[i]] <- hamby252.cmps
  
}

################################
# save data as csv
data_path <- "~/Research/CMPSpaper/CMPSpaper_writing/data/"

for(i in 1:N){
  write.csv(
    CMPS_hamby252_results$cmps.table[[i]] %>% 
      select(-cmps.table, -cmps.table.m) %>% as.data.frame(),
    file = paste(data_path, CMPS_hamby252_results$filename[[i]], ".csv", sep = ""))
}



#################################################
# generate plots
for (i in 1:3) {
  hamby252.cmps <- CMPS_hamby252_results$cmps.table[[i]]
  
  hamby252.plot.list <- list()
  titlee <- CMPS_hamby252_results$titlee[[i]]
  
  hamby252.plot.list[[1]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.max.m,
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
    font("x.text", size = 6)
  
  hamby252.plot.list[[2]] <- hamby252.cmps %>% ggplot() +
    geom_histogram(aes(x = cmps.maxbar.m,
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
    theme(panel.grid.minor = element_blank())
  
  plot <- ggarrange(plotlist = hamby252.plot.list,
                    nrow = 1,
                    ncol = 2,
                    common.legend = TRUE, legend = "bottom")
  plot <- annotate_figure(plot, 
                          top = text_grob(com.title252))
  CMPS_hamby252_results$plot[[i]] <- plot
}

CMPS_hamby252_results$plot[[1]]
CMPS_hamby252_results$plot[[2]]
CMPS_hamby252_results$plot[[3]]

####
# check 281 for second setup
tt1 <- CMPS_hamby252_results$cmps.table[[2]]$cmps.max.m
tt2 <- hamby252_results$cmps.table[[4]]$cmps.max.m
ck_idx <- tt1 != tt2
(1:595)[ck_idx]
tt1 - tt2
tt1[ck_idx]


### phase
# 1 2 3 4 5 6
# 2 3 4 5 6 1
# 3 4 5 6 1 2


(b252$sigs75[[1]] %>% ggplot() +
  geom_line(aes(x=x, y=sig))) / 
(b252$sigs25[[1]] %>% ggplot() +
  geom_line(aes(x=x, y=sig))) / 
(b252$sigs15[[1]] %>% ggplot() +
   geom_line(aes(x=x, y=sig)))

(b252$sigs75[[1]] %>% ggplot() +
    geom_line(aes(x=x, y=value-raw_sig))) / 
  (b252$sigs25[[1]] %>% ggplot() +
     geom_line(aes(x=x, y=value-raw_sig))) / 
  (b252$sigs15[[1]] %>% ggplot() +
     geom_line(aes(x=x, y=value-raw_sig)))




