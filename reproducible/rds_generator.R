###############################
# Note: Please set the working 
# directory so that "./data-csv/"
# is available
###############################

###############################
# for the random forest model
rf.results <- list()
rf.results[[1]] <- read.csv(file = "./data-csv/hamby252/hamby252_rf_results.csv")
rf.results[[2]] <- read.csv(file = "./data-csv/hamby44/hamby44_rf_results.csv")

###############################
# hamby 252 segments
N <- 8 # the number of different sets of parameters
CMPS_hamby252_results <- list() # a container for everything

#### setup span1
CMPS_hamby252_results$span1 <- as.list(rep(0.25, N))

#### setup signature name
CMPS_hamby252_results$signame <- as.list(rep("sigs25", N))
# list("sigs75", "sigs25", "sigs25_1062", "sigs15")

#### setup npeaks.set
CMPS_hamby252_results$npeaks.set <- lapply(1:N, function(t) c(5,3,1))

#### setup seg_length
CMPS_hamby252_results$seg_length <- as.list(c(25, 50, 75, 100, 125, 150, 175, 200))

#### setup Tx
CMPS_hamby252_results$Tx <- as.list(rep(25, N))

#### setup outlength
CMPS_hamby252_results$outlength <- vector(mode = "list", length = N)

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
      "hamby252_segment",
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


for(i in 1:N){
  CMPS_hamby252_results$cmps.table[[i]] <- 
    read.csv(file = paste("./data-csv/hamby252/", 
                          CMPS_hamby252_results$filename[[i]], 
                          ".csv",
                          sep = ""))
}

CMPS_hamby252_results_seg <- CMPS_hamby252_results

###############################
# hamby 252 npeaks
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


for(i in 1:N){
  CMPS_hamby252_results$cmps.table[[i]] <- 
    read.csv(file = paste("./data-csv/hamby252/", 
                          CMPS_hamby252_results$filename[[i]], 
                          ".csv",
                          sep = ""))
}

CMPS_hamby252_results_npeak <- CMPS_hamby252_results

###############################
# hamby 44 segments
N <- 8 # the number of different sets of parameters
CMPS_hamby44_results <- list() # a container for everything

#### setup span1
CMPS_hamby44_results$span1 <- as.list(rep(0.25, N))

#### setup signature name
CMPS_hamby44_results$signame <-
  as.list(rep("sigs25_531", N))

#### setup neapks.set
CMPS_hamby44_results$npeaks.set <-
  lapply(1:N, function(t) c(5,3,1))

#### setup seg_length
CMPS_hamby44_results$seg_length <- as.list(c(30, 45, 61, 90, 122, 150, 180, 210))

#### setup outlength
CMPS_hamby44_results$outlength <- vector(mode = "list", length = N)

#### setup Tx
CMPS_hamby44_results$Tx <- as.list(rep(30, N))

#### create a plot title and a file name for each set of parameters
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
      "hamby44_segment",
      CMPS_hamby44_results$span1[[i]]*100,
      paste(CMPS_hamby44_results$npeaks.set[[i]], collapse = "-"),
      CMPS_hamby44_results$seg_length[[i]],
      CMPS_hamby44_results$Tx[[i]],
      sep = "_"
    )
}

#### container for the cmps results and corresponding plots
CMPS_hamby44_results$cmps.table <- list()
CMPS_hamby44_results$plot <- list()


for(i in 1:N){
  CMPS_hamby44_results$cmps.table[[i]] <- 
    read.csv(file = paste("./data-csv/hamby44/", 
                          CMPS_hamby44_results$filename[[i]], 
                          ".csv",
                          sep = ""))
}

CMPS_hamby44_results_seg <- CMPS_hamby44_results

###############################
# hamby 44 segments
N <- 8 # the number of different sets of parameters
CMPS_hamby44_results <- list() # a container for everything

#### setup span1
CMPS_hamby44_results$span1 <- as.list(rep(0.25, N))

#### setup signature name
CMPS_hamby44_results$signame <-
  as.list(rep("sigs25_531", N))

#### setup neapks.set
CMPS_hamby44_results$npeaks.set <-
  list(c(5),
       c(5,2),
       c(5,3,1),
       c(10,6,2),
       c(6,4,2,1),
       c(10,7,4,2),
       c(10,6,4,2,1),
       c(10,8,6,4,2,1))

#### setup seg_length
CMPS_hamby44_results$seg_length <- as.list(rep(61, N))

#### setup outlength
CMPS_hamby44_results$outlength <- vector(mode = "list", length = N)
CMPS_hamby44_results$outlength[[7]] <- c(61, 91, 121, 151, 181)
CMPS_hamby44_results$outlength[[8]] <- c(61, 91, 121, 151, 181, 211)

#### setup Tx
CMPS_hamby44_results$Tx <- as.list(rep(30, N))

#### create a plot title and a file name for each set of parameters
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
      "hamby44_npeak",
      CMPS_hamby44_results$span1[[i]]*100,
      paste(CMPS_hamby44_results$npeaks.set[[i]], collapse = "-"),
      CMPS_hamby44_results$seg_length[[i]],
      CMPS_hamby44_results$Tx[[i]],
      sep = "_"
    )
}

#### container for the cmps results and corresponding plots
CMPS_hamby44_results$cmps.table <- list()
CMPS_hamby44_results$plot <- list()

for(i in 1:N){
  CMPS_hamby44_results$cmps.table[[i]] <- 
    read.csv(file = paste("./data-csv/hamby44/", 
                          CMPS_hamby44_results$filename[[i]], 
                          ".csv",
                          sep = ""))
}

CMPS_hamby44_results_npeak <- CMPS_hamby44_results

all_results_container <- list(
  h252_npeak = CMPS_hamby252_results_npeak,
  h252_seg = CMPS_hamby252_results_seg,
  h44_npeak = CMPS_hamby44_results_npeak,
  h44_seg = CMPS_hamby44_results_seg,
  rf = rf.results
)

saveRDS(all_results_container, "./CMPSpaper_results.rds")

