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
compute_avgscore_denoise <- function (land1, land2, score, addNA = FALSE) 
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
  scores <- sapply(scores.list, mean, na.rm=TRUE)
  max.phase <- which.max(scores)
  result.match <- max(scores, na.rm = TRUE)
  result.nmatch <- 
    scores.list[-max.phase] %>% unlist() %>% mean(na.rm=TRUE)
  result.match - result.nmatch
}