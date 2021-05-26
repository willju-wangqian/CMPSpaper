p <- 1:595

cb.idx <- 442

tmp.lands <-
  c(
    b44 %>% filter(bullet_id %in% b.cb[, cb.idx][1]) %>% pull(scan_id),
    b44 %>% filter(bullet_id %in% b.cb[, cb.idx][2]) %>% pull(scan_id)
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
      sig1 <- b44$sigs_main[b44$scan_id == xx][[1]]$sig
      sig2 <- b44$sigs_main[b44$scan_id == yy][[1]]$sig
      
      extract_feature_cmps(
        sig1,
        sig2,
        npeaks.set = CMPS_hamby44_results$npeaks.set[[i]],
        seg_length = CMPS_hamby44_results$seg_length[[i]],
        Tx = CMPS_hamby44_results$Tx[[i]],
        include = NULL
      )
    }
  )

for(j in 1:36) {
  print(j) 
  xx <- tmp.comp$land1[j]
  yy <- tmp.comp$land2[j]
  
  sig1 <- b44$sigs_main[b44$scan_id == xx][[1]]$sig
  sig2 <- b44$sigs_main[b44$scan_id == yy][[1]]$sig
  
  extract_feature_cmps(
    sig1,
    sig2,
    npeaks.set = CMPS_hamby44_results$npeaks.set[[i]],
    seg_length = CMPS_hamby44_results$seg_length[[i]],
    Tx = CMPS_hamby44_results$Tx[[i]],
    include = NULL
  )
}

j <- 6

segs <- get_segs(sig1, len = 61)
cmps.result.test <- extract_feature_cmps(sig1, sig2, seg_length = 61, Tx = 30, include = NULL)

seg <- segs$segs[[1]]
comp <- sig2

min.overlap <- min(
  length(seg[!is.na(seg)])*0.9, 
  round(length(comp[!is.na(comp)])*0.1)
)
x <- comp
y <- seg
ccr2 <- get_ccf5(x, y, min.overlap = 62)
ccr2$ccf %>% tail(100)
length(ccr2$ccf)

ccr <- get_ccf4(x, y, min.overlap = 62)
ccr$ccf %>% tail(100)
length(ccr$ccf)









ccr <- get_ccf4(comp, seg, min.overlap = 1)
ccr <- get_ccf5(comp, seg, min.overlap = 0)

CMPS::get_ccf5(comp, seg, min.overlap = 1) -> tt.result

na_trim_c <- CMPS:::na_trim_c
compute_cross_corr_c <- CMPS:::compute_cross_corr_c

min.overlap <- 0
y <- c(y, 1, 1)



find_maxs <- zoo::rollapply(ccr$ccf, 3, function(x) max(x) == x[2],
                       fill = list(NA, NA, NA))
peaks <- which(find_maxs)

peaks <- NULL
peaks <- integer(0)
npeaks <- 5
od <- order(ccr$ccf[peaks], decreasing = TRUE)[1:npeaks]
# adjust the position
adj_pos <- ccr$lag - segs$index[[3]][1] + 1

peaks.heights <- ccr$ccf[peaks][od]
peaks.pos <- adj_pos[peaks][od]

peaks.pos <- c(peaks.pos, 1)
ck.tmp <- abs(peaks.pos - 5) <= 31 
peaks.pos[ck.tmp]



x.na.count <- .Call(na_trim_c, x)
# .Call("_NA_TRIM", x)
# y.na.count <- na_trim_cmps(y) #.Call("_NA_TRIM", y)
y.na.count <- .Call(na_trim_c, y)

x.narm <- x[(x.na.count[1] + 1) : (length(x) - x.na.count[2])]
y.narm <- y[(y.na.count[1] + 1) : (length(y) - y.na.count[2])]
# nxx <- length(x.narm)
nyy <- length(y.narm)

xx <- c(rep(NA, nyy - min.overlap), x.narm, rep(NA, nyy - min.overlap))
# yy <- c(y, rep(NA, length(xx) - ny))

# make sure that y has no NA in the front or in the end
# make sure that xx has no extra NA values; all NA are needed for shifting
# make sure that ny >= min.overlap
# assume x and y have been na.trim-ed
# cors <- compute_cross_corr(xx, y.narm, as.numeric(min.overlap))
cors <- .Call(compute_cross_corr_c, xx, y.narm, as.numeric(min.overlap))
cors <- c(rep(NA, x.na.count[1] + y.na.count[2]), 
          cors, 
          rep(NA, x.na.count[2] + y.na.count[1]))
length(cors)


for(cb.idx in 495:101) {
  print(cb.idx)
  
  tmp.lands <-
    c(
      b44 %>% filter(bullet_id %in% b.cb[, cb.idx][1]) %>% pull(scan_id),
      b44 %>% filter(bullet_id %in% b.cb[, cb.idx][2]) %>% pull(scan_id)
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
        sig1 <- b44$sigs_main[b44$scan_id == xx][[1]]$sig
        sig2 <- b44$sigs_main[b44$scan_id == yy][[1]]$sig
        
        extract_feature_cmps(
          sig1,
          sig2,
          npeaks.set = CMPS_hamby44_results$npeaks.set[[i]],
          seg_length = CMPS_hamby44_results$seg_length[[i]],
          Tx = CMPS_hamby44_results$Tx[[i]],
          include = "nseg"
        )
      }
    )
  
  
  
  # cat(cb.idx, "- 3; ")
  
  # tmp.comp <- tmp.comp %>%
  #   mutate(
  #     cmps_score = sapply(tmp.comp$cmps, function(x)
  #       x$CMPS.score),
  #     cmps_nseg = sapply(tmp.comp$cmps, function(x)
  #       x$nseg)
  #   )
  # 
  # # cat(cb.idx, "- 4;\n")
  # 
  # cmps.table <-
  #   tmp.comp %>% select(land1, land2, cmps_score, cmps_nseg)
  # # cp2 <- tmp.comp %>% select(land1, land2, cmps_score, cmps_nseg)
  # 
  # # cat(" end of ", cb.idx, "\n", "############# \n")
  # 
  # tibble(
  #   bullet1 = b.cb[, cb.idx][1],
  #   bullet2 = b.cb[, cb.idx][2],
  #   cmps.table = list(cmps.table)
  # )
}



tmp.44.list <- parLapply(cl, p, function(cb.idx) {
  # cat(" ############# \n", "start of ", cb.idx, "\n")
  
  tmp.lands <-
    c(
      b44 %>% filter(bullet_id %in% b.cb[, cb.idx][1]) %>% pull(scan_id),
      b44 %>% filter(bullet_id %in% b.cb[, cb.idx][2]) %>% pull(scan_id)
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
        sig1 <- b44$sigs_main[b44$scan_id == xx][[1]]$sig
        sig2 <- b44$sigs_main[b44$scan_id == yy][[1]]$sig
        
        extract_feature_cmps(
          sig1,
          sig2,
          npeaks.set = CMPS_hamby44_results$npeaks.set[[i]],
          seg_length = CMPS_hamby44_results$seg_length[[i]],
          Tx = CMPS_hamby44_results$Tx[[i]],
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
})