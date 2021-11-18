We thank you for the time and efforts you have invested in the review.

## Code

> library(bulletxtrctr) is not available on CRAN, and no instructions are given for compiling. It is needed for bullet_to_land_predict() later in the code. Technically, the requirement that the package is on CRAN also means that all the dependent packages are also on CRAN. Is there a reason for bulletxtrctr not to be on CRAN? If so, then your example code for cmpsR will need to be independent of bulletxtrctr.

something

> source("../code/func_collection.R") points to a non-existing directory
>
> This function compute_var_ratio_anova() isn’t found.

Thank you for pointing this out. This was an uncautious mistake and has been fixed. `../code/func_collection.R` is not needed for the paper. `compute_var_ratio_anova()` should actually be `cmpsR::compute_ss_ratio()` and is part of the `cmpsR` package.

## Paper

> If the paper was written in Rmarkdown please submit the .Rmd file too. My recommendation is for you to follow the format provided in the (new) rjtools package (https://rjournal.github.io/rjtools/), when writing with Rmarkdown. (You should probably remove the extra lines referring to knitr options from the .R file after using curl to extract code from the Rmd.) (More details on this format are provided on the dev version of the R Journal site, which will become the primary web site from the next issue: https://journal.r-project.org/dev/submissions.html .)

Yes, the paper was written in Rmarkdown, and the `.Rmd` file has been included in the re-submission. 

The extra lines referring to knitr options from the .R file are removed. 

To address the requirements described on the dev version of the R Journal site, all the R scripts in the `reproducible` folder of the original submission has been zipped into `supplementary-files.zip`, and the data contained in the `reproducible` folder has been made available on an alternative site with the link provided in the paper under the section "Supplementary materials" to reduce the size of the submission file.

> Please don’t have URLs randomly in the text: “Our implementation is made available as part of the R package cmpsR on GitHub at https://github.com/willju-wangqian/ cmpsR.” You should be referring to the CRAN version, and development version on github, and these would be listed in a section titled “Supplementary materials” after the Conclusions section.

> Sentences like “(We have asked the authors of Chen et al. (2019), and we are still waiting for the replies.)” should be omitted from the paper.

This sentence has been removed from the paper. 

> The paper is nicely written and the plots are easy to read.

Thank you for
