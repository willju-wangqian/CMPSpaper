library(assertthat)

compute_var_ratio <- function(score, label) {
  score.split <- split(score, label)
  
  group.mean <- sapply(score.split, mean, na.rm=TRUE)
  
  total.mean <- mean(score, na.rm = TRUE)
  
  var.between <- sum((group.mean - total.mean)^2)
  var.within <- sum(sapply(score.split, var, na.rm = TRUE))
  
  var.ratio <- var.between / var.within
  return(var.ratio)
}

compute_var_ratio_anova <- function(score, label, MS=TRUE) {
  score.split <- split(score, label)
  
  group.mean <- sapply(score.split, mean, na.rm=TRUE)
  group.count <- sapply(score.split, function(tt) sum(!is.na(tt)) )
  total.mean <- mean(score, na.rm = TRUE)
  
  k <- length(unique(label))
  n <- length(score)
  
  var.between <- sum((group.mean - total.mean)^2 * group.count)
  var.within <- sum(sapply(score.split, function(tt) {
    sum((tt - mean(tt, na.rm = TRUE))^2)
  })) 
  
  var.ratio <- var.between / var.within
  if(MS) {
    var.ratio <- var.ratio * (n-k) / (k-1)
  }
  return(var.ratio)
}

get_all_phases <- function(land1, land2, score, addNA = FALSE)
{
  maxland <- max(land1, land2)
  fullframe <- data.frame(expand.grid(land1 = 1:maxland, land2 = 1:maxland))
  bcompare <- data.frame(land1, land2, score)
  fullframe <- fullframe %>% left_join(bcompare, by = c("land1", 
                                                        "land2"))
  fullframe <- fullframe %>% mutate(land1 = factor(land1, levels = 1:maxland), 
                                    land2 = factor(land2, levels = 1:maxland))
  matrix <- xtabs(score ~ land1 + land2, data = fullframe, 
                  addNA = addNA)/xtabs(~land1 + land2, data = fullframe, 
                                       addNA = addNA)
  matrix <- cbind(matrix, matrix)
  scores.list <- 1:maxland %>% lapply(FUN = function(i) {
    if (i == 1) {
      diag(matrix)
    }
    else {
      i <- i - 1
      diag(matrix[, -(1:i)])
    }
  })
  
  return(scores.list)
}

## the function re-written from compute_average_scores
compute_avgscore_denoise <- function (land1, land2, score, FUNC = mean, addNA = FALSE, na.rm = TRUE) 
{
  
  scores.list <- get_all_phases(land1, land2, score, addNA)
  
  result <- compute_diff_phase(scores.list, FUNC, na.rm)
  
  return(result)
}

compute_diff_phase <- function(scores.list, FUNC = mean, na.rm = TRUE, both = FALSE)
{
  scores <- sapply(scores.list, FUNC, na.rm=na.rm)
  max.phase <- which.max(scores)
  result.match <- max(scores, na.rm=na.rm)
  result.nmatch <- scores.list[-max.phase] %>% unlist() %>% FUNC(na.rm=na.rm)
  if(both) {
    return(c(result.match, result.nmatch))
  } else {
    return(result.match - result.nmatch)
  }
}

cmps_metrics_helper <- function(cmps.table) {
  
  landidx1 <- cmps.table$landidx1
  landidx2 <- cmps.table$landidx2
  cmps_score <- cmps.table$cmps_score
  cmps_score_scaled <- cmps.table$cmps_score_scaled
  result <- data.frame(
    cmps.diff = compute_avgscore_denoise(landidx1, landidx2, cmps_score, addNA = TRUE),
    cmps.diff.med = compute_avgscore_denoise(landidx1, landidx2, 
                                              cmps_score, FUNC = median, addNA = TRUE),
    cmps.max = max(cmps_score, na.rm = TRUE),
    cmps.maxbar = compute_average_scores(landidx1, landidx2, cmps_score) %>% max(na.rm = TRUE),
    cmps.diff_scaled = compute_avgscore_denoise(landidx1, landidx2, cmps_score_scaled, addNA = TRUE),
    cmps.diff.med_scaled = compute_avgscore_denoise(landidx1, landidx2, 
                                                    cmps_score_scaled, FUNC = median, addNA = TRUE),
    cmps.max_scaled = max(cmps_score_scaled, na.rm = TRUE),
    cmps.maxbar_scaled = compute_average_scores(landidx1, landidx2, cmps_score_scaled) %>% 
      max(na.rm = TRUE)
  )
  # result$cmps.diff <- compute_avgscore_denoise(landidx1, landidx2, cmps_score, addNA = TRUE)
  # result$cmps.diff.med <- compute_avgscore_denoise(landidx1, landidx2, 
  #                                                  cmps_score, FUNC = median, addNA = TRUE)
  # result$cmps.max <- max(cmps_score, na.rm = TRUE)
  # result$cmps.maxbar <- compute_average_scores(landidx1, landidx2, cmps_score) %>% max(na.rm = TRUE)
  # 
  # result$cmps.diff_scaled <- compute_avgscore_denoise(landidx1, landidx2, cmps_score_scaled, addNA = TRUE)
  # result$cmps.diff.med_scaled <- compute_avgscore_denoise(landidx1, landidx2, 
  #                                                  cmps_score_scaled, FUNC = median, addNA = TRUE)
  # result$cmps.max_scaled <- max(cmps_score_scaled, na.rm = TRUE)
  # result$cmps.maxbar_scaled <- 
  #   compute_average_scores(landidx1, landidx2, cmps_score_scaled) %>% max(na.rm = TRUE)
  return(result)
}

compute_score_metrics <- function(land1, land2, score, 
                                  addNA = TRUE, na.rm = TRUE, include = NULL, out_names = NULL) {
  
  assert_that(
    is.numeric(land1), is.numeric(land2), is.numeric(score)
  )
  
  scores.list <- get_all_phases(land1, land2, score, addNA = addNA)
  tmp.diff <- compute_diff_phase(scores.list, mean, na.rm, both = TRUE)
  
  result <- data.frame(
    diff = tmp.diff[1] - tmp.diff[2],
    diff.med = compute_diff_phase(scores.list, median, na.rm = na.rm),
    max = max(score, na.rm = na.rm),
    maxbar = tmp.diff[1]
  )
  
  
  if(is.null(include)) {
    if(!is.null(out_names) && length(out_names) == ncol(result)) {
      names(result) <- out_names
    } else if (!is.null(out_names)) {
      warning("Warning: Fail to change the variable names in the result.")
    }
    return(result)
  } else {
    mm <- match.arg(include, c("diff", "diff.med", "max", "maxbar"), several.ok = TRUE)
    
    if(!is.null(out_names) && length(out_names) == length(mm)) {
      names(result[mm]) <- out_names
      return(result[out_names])
    } else if (!is.null(out_names)) {
      warning("Warning: Fail to change the variable names in the result.")
    }
    
    return(result[mm])
    
  }
}


score_metrics_helper <- function(cmps.table, score = "cmps_score", 
                                 scaled = NULL, addNA = TRUE, na.rm = TRUE) {
  
  landidx1 <- cmps.table$landidx1
  landidx2 <- cmps.table$landidx2
  score <- cmps.table[[score]]
  result <- data.frame(
    diff = compute_avgscore_denoise(landidx1, landidx2, score, addNA = addNA),
    diff.med = compute_avgscore_denoise(landidx1, landidx2, 
                                             score, FUNC = median, addNA = addNA),
    max = max(score, na.rm = na.rm),
    maxbar = compute_average_scores(landidx1, landidx2, score, addNA = addNA) %>% max(na.rm = na.rm)
  )
  
  if(!is.null(scaled)) {
    score_scaled <- cmps.table[[scaled]]
    result <- result %>% mutate(
      diff_scaled = compute_avgscore_denoise(landidx1, landidx2, score_scaled, addNA = addNA),
      diff.med_scaled = compute_avgscore_denoise(landidx1, landidx2, 
                                                      score_scaled, FUNC = median, addNA = addNA),
      max_scaled = max(score_scaled, na.rm = na.rm),
      maxbar_scaled = compute_average_scores(landidx1, landidx2, score_scaled, addNA = addNA) %>% 
        max(na.rm = na.rm)
    )
  }

  return(result)
}

metric_plot_helper <- function(cmps.metric, metric, scaled = FALSE, SSratio = TRUE, ...) {
  # scaled <- str_detect(metric, "_scaled")
  dots <- list(...)
  if(scaled) {
    if(is.null(dots$breaks)) {
      dots$breaks <- seq(0,1,0.05)
    }
    
    if(is.null(dots$binwidth)) {
      dots$binwidth <- 0.04
    }
    
  } else {
    if(is.null(dots$breaks)) {
      dots$breaks <- seq(0,24,1)
    }
    
    if(is.null(dots$binwidth)) {
      dots$binwidth <- 1
    }
  }
  
  p <- cmps.metric %>% ggplot() +
    geom_histogram(aes(x = .data[[metric]],
                       fill = as.factor(.data$type_truth)), binwidth = dots$binwidth) +
    labs(
      x = metric,
      fill = "Comparison Type",
      subtitle = dots$subtitle
    ) +
    scale_x_continuous(breaks = dots$breaks) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(values=c("darkorange", "darkgrey"))
  
  ss.ratio <- compute_var_ratio_anova(cmps.metric[[metric]], cmps.metric$type_truth, MS = FALSE)
  if(SSratio) {
    p <- p + 
      annotate(geom = "label", x = Inf, y = Inf, 
               label = paste("SS Ratio:", round(ss.ratio, 2)), 
               fill = "white", hjust = 1, vjust = 1)
  }
  
  return(p)
}


###################
# update sig_align
sig_align <- function (sig1, sig2) 
{
  assertthat::assert_that(is.numeric(sig1), is.numeric(sig2))
  sig1 <- CMPS::cmps_na_trim(sig1)
  sig2 <- CMPS::cmps_na_trim(sig2)
  
  n1 <- length(sig1)
  n2 <- length(sig2)
  if (n1 > n2) {
    x <- sig1
    y <- sig2
  }
  else {
    x <- sig2
    y <- sig1
  }
  cors <- CMPS::get_ccf4(x, y, round(0.75 * min(length(sig1), length(sig2))))
  lag <- cors$lag[which.max(cors$ccf)]
  if (lag < 0) {
    x <- c(rep(NA, abs(lag)), x)
  }
  if (lag > 0) {
    y <- c(rep(NA, lag), y)
  }
  delta <- length(x) - length(y)
  if (delta < 0) 
    x <- c(x, rep(NA, abs(delta)))
  if (delta > 0) 
    y <- c(y, rep(NA, delta))
  if (n1 > n2) {
    dframe0 <- data.frame(x = 1:length(x), sig1 = x, sig2 = y)
  }
  else {
    dframe0 <- data.frame(x = 1:length(x), sig1 = y, sig2 = x)
  }
  maxcor <- max(cors$ccf, na.rm = TRUE)
  list(ccf = maxcor, lag = lag, lands = dframe0)
}


get_ccf2 <- function (x, y, min.overlap = round(0.1 * max(length(x), length(y)))) 
{
  x <- as.vector(unlist(x))
  y <- as.vector(unlist(y))
  nx <- length(x)
  ny <- length(y)
  assert_that(is.numeric(x), is.numeric(y))
  assert_that(nx > 0, ny > 0, nx >= ny)
  xx <- c(rep(NA, ny - min.overlap), x, rep(NA, ny - min.overlap))
  yy <- c(y, rep(NA, length(xx) - ny))
  lag.max <- length(yy) - length(y)
  lags <- 0:lag.max
  cors <- sapply(lags, function(lag) {
    cor(xx, lag(yy, lag), use = "pairwise.complete")
  })
  ns <- sapply(lags, function(lag) {
    dim(na.omit(cbind(xx, lag(yy, lag))))[1]
  })
  cors[ns < min.overlap] <- NA
  lag <- lags - (ny - min.overlap)
  return(list(lag = lag, ccf = cors))
}
  
# aligned <- tmp.comp$aligned[[1]]$lands
my_extract_feature_lag <- function(aligned)
{
  assert_that(dim(aligned)[2] > 2, msg = "aligned must have at least 3 columns")
  for (i in 2:dim(aligned)[2]) {
    assert_that(is.numeric(aligned[, i]), msg = sprintf("Column %d (%s) is not numeric", 
                                                        i, names(aligned)[i]))
  }
  lags <- sapply(aligned[, -1], function(x) {
    if (!is.na(x[1])) 
      return(0)
    diffs <- diff(is.na(x))
    which(diffs == -1)
  })
  if (length(lags) == 2) 
    return(diff(lags))
  lags
}

my_extract_feature_all <- function (aligned, striae, resolution, tmpfile = NULL, ...) 
{
  feature <- value <- NULL
  assert_that(!is.null(aligned), !is.null(striae), msg = "aligned and striae must not be NULL")
  features <- apropos("^extract_feature_") # change: added a start anchor
  features <- features[!str_detect(features, "extract_feature_cmps")] # change: exclude _cmps from the search list
  dots <- list(...)
  values <- features %>% purrr::map_dbl(.f = function(f) {
    fun <- getFromNamespace(f, asNamespace("bulletxtrctr"))
    fun_args <- names(formals(fun))
    matching_args <- dots[names(dots) %in% fun_args]
    if ("aligned" %in% fun_args) {
      matching_args$aligned <- aligned$lands
    }
    if ("striae" %in% fun_args) {
      matching_args$striae <- striae$lines
    }
    if ("resolution" %in% fun_args) {
      matching_args$resolution <- resolution
    }
    res <- do.call(fun, matching_args)
    res
  })
  dframe <- data.frame(feature = gsub("extract_feature_", "", 
                                      features), value = values) %>% spread(feature, value)
  if (!is.null(tmpfile)) {
    if (file.exists(tmpfile)) {
      write.table(dframe, file = tmpfile, sep = ",", append = TRUE, 
                  col.names = FALSE, row.names = FALSE)
    }
    else {
      write.table(dframe, file = tmpfile, sep = ",", append = FALSE, 
                  col.names = TRUE, row.names = FALSE)
    }
  }
  dframe
}

