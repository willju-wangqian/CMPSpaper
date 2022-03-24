The `supplementary-files.zip` contains R scripts for reproducing all aspects of the paper "An Open-Source Implementation of the CMPS Algorithm for Assessing Similarity of Bullets"

**Data not included in "supplementary-files.zip"**

The processed data of Hamby set 252 and Hamby set 44 are stored in two `.rds` files:

-   [`BulletSignatures252.rds`](https://github.com/willju-wangqian/CMPSpaper/blob/main/reproducible/bullet_signatures_etc/BulletSignatures252.rds)
-   [`BulletSignatures44.rds`](https://github.com/willju-wangqian/CMPSpaper/blob/main/reproducible/bullet_signatures_etc/BulletSignatures44.rds)

And the ground truth of Hamby set 252 is provided in [`StudyInfo.xlsx`](https://github.com/willju-wangqian/CMPSpaper/blob/main/reproducible/bullet_signatures_etc/StudyInfo.xlsx)

Due to the size of the file, these processed data are not included in the "supplementary-files.zip", but the links are provided for download.

In order to reproduce the results in the paper, please download "BulletSignatures252.rds", "BulletSignatures44.rds", and "StudyInfo.xlsx" and save them in a folder named "bullet_signatures_etc". And place the folder "bullet_signatures_etc" into the folder of all the R scripts of "supplementary-files.zip".

Please check out the folder structure of [the `reproducible` folder](https://github.com/willju-wangqian/CMPSpaper/tree/main/reproducible) for reference.

**File description**

-   `hamby*_result_generator.R`: these R scripts take the processed data of Hamby set 252 and Hamby set 44, generate preliminary results used in the paper, and save the results in the `.csv` format in the folder `data-csv`.
-   `rds_generator.R`: this R script takes the generated `.csv` files and produce `CMPSpaper_results.rds`. `CMPSpaper_results.rds` is identical to `data/CMPSpaper_results.rds` and is used generate figures and other results presented in the paper.
-   `func_collection.R`: this R script contains some helper functions used by other R scripts.

**Data included in "supplementary-files.zip"**

-   `data-csv`: this folder contains `.csv` files we generated using R scripts `hamby*_result_generator.R` 
-   `CMPSpaper_results.rds`: this is the `.rds` file we generated using R script `rds_generator.R`

These data can be used as reference or example results of the reproducible codes

-   `csafe_rf2.rds`: this `.rds` file contains the random forest model used in this paper

