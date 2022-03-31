# run all scripts
source("./hamby44_seg_result_generator.R")
source("./hamby44_npeak_result_generator.R")
source("./hamby44_rf_result_generator.R")

source("./hamby252_seg_result_generator.R")
source("./hamby252_npeak_result_generator.R")
source("./hamby252_rf_result_generator.R")

# test if we get the same result
new_files_name <- list.files("./data-csv", recursive = T, full.names = T)
ogn_files_name <- list.files("./data-csv-Copy", recursive = T, full.names = T)

for (i in 1:length(new_files_name)) {
  tt1 <- read.csv(new_files_name[i])
  tt2 <- read.csv(ogn_files_name[i])
  
  if(!identical(tt1,tt2)) print(i)
}

# 5, 22
i <- 22
new_files_name[i]
ogn_files_name[i]
