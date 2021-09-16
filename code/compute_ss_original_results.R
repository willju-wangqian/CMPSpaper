library(magrittr, include.only = c("%<>%", "%$%"))

result.container <- list()

result.container$npeak <- c("n5", "n531", "n5", "n531")
result.container$metrics <- c("cmps.max", "cmps.max", "cmps.maxbar", "cmps.maxbar")

dd.list <- list()

dd.list[[1]] <- tibble(
  score = c(9:13, 16:24),
  freq = c(132, 282, 115, 19, 1, 1, 3,2,4,8,12,7,7,2),
  type = c(rep("KNM", 5), rep("KM", 9))
)

dd.list[[2]] <- tibble(
  score = c(1:5, 14:23),
  freq = c(7,343,168,24,7,1,1,2,3,6,7,8,8,7,3),
  type = c(rep("KNM", 5), rep("KM", 10))
)

dd.list[[3]] <- tibble(
  score = c(6:9, 13:19),
  freq = c(2,210,307,30,6,3,9,9,11,7,1),
  type = c(rep("KNM", 4), rep("KM", 7))
)

dd.list[[4]] <- tibble(
  score = c(0:2, 8:17),
  freq = c(36, 508, 5, 1,3,4,2,11,3,9,8,4,1),
  type = c(rep("KNM", 3), rep("KM", 10))
)

result.container$dd.list <- dd.list

result.container <- result.container %>% as_tibble()

result.container %<>% as_tibble()

result.container %<>% mutate(
  dd.table = purrr::map(dd.list, function(dd) {
    tibble(score = rep(dd$score, dd$freq),
           type_truth = rep(dd$type, dd$freq))
  }),
  ss.ratio = purrr::map_dbl(dd.table, function(dd) {
    compute_var_ratio_anova(dd$score, dd$type_truth, MS = FALSE)
  })
) 










