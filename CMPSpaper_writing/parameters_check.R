N <- 3
CMPS_hamby44_results <- list()
CMPS_hamby44_results$span1 <- list(0.25, 0.25, 0.75)
CMPS_hamby44_results$signame <-
  list("sigs25_531", "sigs25_1061", "sigs75_531")
CMPS_hamby44_results$npeaks.set <-
  list(c(5, 3, 1),
       c(10, 6, 1),
       c(5, 3, 1))
CMPS_hamby44_results$seg_length <- as.list(rep(61, 3))
CMPS_hamby44_results$Tx <- as.list(rep(30, 3))
CMPS_hamby44_results$titlee <- list()
CMPS_hamby44_results$filename <- list()
for (i in 1:N) {
  CMPS_hamby44_results$titlee[[i]] <-
    paste0(
      "npeaks.set=c(",
      paste(CMPS_hamby44_results$npeaks.set[[i]], collapse = ","),
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
      "hamby44",
      CMPS_hamby44_results$span1[[i]]*100,
      paste(CMPS_hamby44_results$npeaks.set[[i]], collapse = "-"),
      CMPS_hamby44_results$seg_length[[i]],
      CMPS_hamby44_results$Tx[[i]],
      sep = "_"
    )
}
CMPS_hamby44_results$cmps.table <- list()
CMPS_hamby44_results$plot <- list()

for(i in 1:N){
  CMPS_hamby44_results$cmps.table[[i]] <- 
    read.csv(file = paste("./CMPSpaper_writing/data/", 
                          CMPS_hamby44_results$filename[[i]], 
                          ".csv",
                          sep = ""))
}

tab1 <- CMPS_hamby44_results$cmps.table[[1]] %>% select(cmps.maxbar, type_truth)
anova(tab1)

library(MASS)

tab1 <- CMPS_hamby44_results$cmps.table[[1]] %>% dplyr::select(cmps.maxbar, type_truth)
tab1$cmps.maxbar.scaled <- (tab1$cmps.maxbar - mean(tab1$cmps.maxbar))/sd(tab1$cmps.maxbar)
model <- lda(type_truth ~ cmps.maxbar.scaled, data = tab1)
model$svd

model.aov <- aov(cmps.maxbar.scaled ~ type_truth, data = tab1)
summary(model.aov)

model1 <- lda(type_truth ~ cmps.maxbar, data = CMPS_hamby44_results$cmps.table[[1]])
model2 <- lda(type_truth ~ cmps.maxbar, data = CMPS_hamby44_results$cmps.table[[2]])
model3 <- lda(type_truth ~ cmps.maxbar, data = CMPS_hamby44_results$cmps.table[[3]])


model1$svd
model2$svd
model3$svd








