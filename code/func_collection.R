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

## the function re-written from compute_average_scores
compute_avgscore_denoise <- function (land1, land2, score, FUNC = mean, addNA = FALSE, na.rm = TRUE) 
{
  # if (!is.numeric(land1)) 
  #   land1 <- readr::parse_number(as.character(land1))
  # if (!is.numeric(land2)) 
  #   land2 <- readr::parse_number(as.character(land2))
  # assert_that(is.numeric(land1), is.numeric(land2), is.numeric(score))
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
  scores <- sapply(scores.list, FUNC, na.rm=na.rm)
  max.phase <- which.max(scores)
  result.match <- max(scores, na.rm=na.rm)
  result.nmatch <- 
    scores.list[-max.phase] %>% unlist() %>% FUNC(na.rm=na.rm)
  result.match - result.nmatch
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

cmps_metric_plot_helper <- function(cmps.metric, metric) {
  scaled <- str_detect(metric, "_scaled")
  if(scaled) {
    p <- cmps.metric %>% ggplot() +
      geom_histogram(aes(x = .data[[metric]],
                         fill = as.factor(.data$type_truth)), binwidth = 0.04) +
      labs(
        x = metric,
        fill = "Comparison Type",
        subtitle = titlee
      ) +
      scale_x_continuous(breaks = seq(0, 1, 0.05)) +
      theme_bw() +
      theme(panel.grid.minor = element_blank()) +
      scale_fill_manual(values=c("darkorange", "darkgrey"))
  } else {
    p <- cmps.metric %>% ggplot() +
      geom_histogram(aes(x = .data[[metric]],
                         fill = as.factor(.data$type_truth)), binwidth = 1) +
      labs(
        x = metric,
        fill = "Comparison Type",
        subtitle = titlee
      ) +
      scale_x_continuous(breaks = seq(0, 24, 1)) +
      theme_bw() +
      theme(panel.grid.minor = element_blank()) +
      scale_fill_manual(values=c("darkorange", "darkgrey"))
  }
  return(p)
}

