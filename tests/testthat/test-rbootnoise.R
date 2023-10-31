#The purpose of these rbootnoise tests is to perform rough comparisons tolerating relatively large deviations from the earlier-acquired reference data. The tests try to catch possible absurd deviations caused by significant technical issues. Without the highly controlled containers (not applicable in the context of cross-platform R CMD Checks) the exact technical reproducibility cannot be established. For example, it is known that even the set.seed() can give varying outcomes depending on the version of R on the same underlying system, an unavoidable technical curiosity accepted by the community. https://stackoverflow.com/questions/47199415/is-set-seed-consistent-over-different-versions-of-r-and-ubuntu 

library(lme4, quietly = TRUE)

data(jsp728, package = "lmeresampler")

# ==============================================================================
context("rbootnoise")
# ==============================================================================

test_that("Compare rbootnoise = 0 to lmeresampler 0.2.2 results before the implementation of the feature",{

  #Reference data creator: Ilmari Tamminen
  
  load("./reference_data_for_rbootnoise_test_A071022.RData")
  
  #Obtained with the following specs:
  
  #> sessionInfo()
  #R version 4.2.1 (2022-06-23 ucrt)
  #Platform: x86_64-w64-mingw32/x64 (64-bit)
  #Running under: Windows 10 x64 (build 19044)
  
  #Matrix products: default
  
  #locale:
  #[1] LC_COLLATE=Finnish_Finland.utf8  LC_CTYPE=Finnish_Finland.utf8    LC_MONETARY=Finnish_Finland.utf8
  #[4] LC_NUMERIC=C                     LC_TIME=Finnish_Finland.utf8    
  
  #attached base packages:
  #[1] stats     graphics  grDevices utils     datasets  methods   base     
  
  #other attached packages:
  #[1] lmeresampler_0.2.2 lme4_1.1-29        Matrix_1.4-1      
  
  #loaded via a namespace (and not attached):
  #[1] Rcpp_1.0.8.3         lubridate_1.8.0      lattice_0.20-45      tidyr_1.2.1          prettyunits_1.1.1   
  #[6] ps_1.7.1             rprojroot_2.0.3      digest_0.6.29        utf8_1.2.2           R6_2.5.1            
  #[11] plyr_1.8.7           HLMdiag_0.5.0        ggplot2_3.3.6        pillar_1.8.1         rlang_1.0.6         
  #[16] curl_4.3.2           rstudioapi_0.13      minqa_1.2.4          callr_3.7.0          nloptr_2.0.3        
  #[21] desc_1.4.1           diagonals_6.4.0      devtools_2.4.3       splines_4.2.1        statmod_1.4.37      
  #[26] stringr_1.4.1        munsell_0.5.0        compiler_4.2.1       janitor_2.1.0        pkgconfig_2.0.3     
  #[31] pkgbuild_1.3.1       mgcv_1.8-40          tidyselect_1.2.0     tibble_3.1.7         fansi_1.0.3         
  #[36] crayon_1.5.1         dplyr_1.0.10         withr_2.5.0          MASS_7.3-57          distributional_0.3.1
  #[41] ggdist_3.2.0         grid_4.2.1           nlme_3.1-158         gtable_0.3.1         lifecycle_1.0.3     
  #[46] magrittr_2.0.3       scales_1.2.1         cli_3.3.0            stringi_1.7.8        cachem_1.0.6        
  #[51] farver_2.1.1         reshape2_1.4.4       fs_1.5.2             remotes_2.4.2        snakecase_0.11.0    
  #[56] ellipsis_0.3.2       generics_0.1.3       vctrs_0.4.1          boot_1.3-28          tools_4.2.1         
  #[61] forcats_0.5.2        rcmdcheck_1.4.0      glue_1.6.2           purrr_0.3.4          processx_3.6.1      
  #[66] pkgload_1.3.0        fastmap_1.1.0        nlmeU_0.70-9         colorspace_2.0-3     xopen_1.0.0         
  #[71] sessioninfo_1.2.2    memoise_2.0.1        usethis_2.1.6
  
  #RStudio 2022.02.0+443 "Prairie Trillium" Release (9f7969398b90468440a501cf065295d9050bb776, 2022-02-16) for Windows
  #Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.8 Chrome/69.0.3497.128 Safari/537.36

  #Obtained with the following code:
  
  #library(lme4)
  #library(lmeresampler)
  #model <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
  #nsim <- 100
  #set.seed(123)
  #A071022ref <- bootstrap(model, .f = fixef, type = "residual", B = nsim)

  #About Ubuntu R-CMD-checks. It was noticed that the results produced with Linux systems have insignificant
  #deviations from the results of other systems. Small relative deviations around 10^-10 were observed between
  #the reference data created using Windows and test data created using Linux. Thus, exactly the same
  #reproduction of the reference data is not expected. Maximum relative deviations of 0.001 between the
  #results created using different operating systems is allowed. The linux-test data (not used in this test_that())
  #on which the above is based on was obtained with the following specs:

  #> sessionInfo()
  #R version 4.2.1 (2022-06-23)
  #Platform: x86_64-pc-linux-gnu (64-bit)
  #Running under: CentOS Linux 7 (Core)

  #Matrix products: default
  #BLAS/LAPACK: /home/opt/easybuild/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib64/libflexiblas.so.3.2

  #locale:
  #[1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
  #[3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
  #[5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
  #[7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
  #[9] LC_ADDRESS=C               LC_TELEPHONE=C
  #[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

  #attached base packages:
  #[1] stats     graphics  grDevices utils     datasets  methods   base

  #other attached packages:
  #[1] lmeresampler_0.2.4 lme4_1.1-29        Matrix_1.4-1

  #loaded via a namespace (and not attached):
  #[1] Rcpp_1.0.8.3         plyr_1.8.7           pillar_1.7.0
  #[4] compiler_4.2.1       nloptr_2.0.3         forcats_0.5.1
  #[7] tools_4.2.1          statmod_1.4.36       boot_1.3-28
  #[10] lubridate_1.8.0      lifecycle_1.0.3      tibble_3.1.7
  #[13] nlme_3.1-158         gtable_0.3.0         lattice_0.20-45
  #[16] mgcv_1.8-40          pkgconfig_2.0.3      HLMdiag_0.5.0
  #[19] rlang_1.0.6          cli_3.6.0            DBI_1.1.3
  #[22] stringr_1.4.0        dplyr_1.0.9          janitor_2.2.0
  #[25] generics_0.1.2       vctrs_0.5.2          diagonals_6.4.0
  #[28] grid_4.2.1           tidyselect_1.1.2     nlmeU_0.70-9
  #[31] snakecase_0.11.0     glue_1.6.2           R6_2.5.1
  #[34] fansi_1.0.3          distributional_0.3.0 minqa_1.2.4
  #[37] tidyr_1.2.0          farver_2.1.0         reshape2_1.4.4
  #[40] ggplot2_3.4.1        purrr_0.3.4          magrittr_2.0.3
  #[43] scales_1.2.0         ellipsis_0.3.2       MASS_7.3-57
  #[46] splines_4.2.1        ggdist_3.2.1         assertthat_0.2.1
  #[49] colorspace_2.0-3     utf8_1.2.2           stringi_1.7.6
  #[52] munsell_0.5.0        crayon_1.5.1
  
  model <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
  nsim <- 100
  set.seed(123)
  boo <- bootstrap(model, .f = fixef, type = "residual", B = nsim)
  
  maxreldev <- (A071022ref[["stats"]][,2:5] - boo[["stats"]][,2:5])/A071022ref[["stats"]][,2:5]
  maxreldev <- max(abs(maxreldev))
  comparison <- (maxreldev < 1)
  expect_true(comparison, info = NULL, label = NULL)

  maxreldev <- (A071022ref[["replicates"]] - boo[["replicates"]])/A071022ref[["replicates"]]
  maxreldev <- max(abs(maxreldev))
  comparison <- (maxreldev < 1)
  expect_true(comparison, info = NULL, label = NULL)
})

test_that("Compare rbootnoise = 0.0001 to the results of the first implementation of the feature",{

  #Reference data creator: Ilmari Tamminen
  
  load("./reference_data_for_rbootnoise_test_B071022.RData")
  
  #Obtained with the following specs:
  
  #> sessionInfo()
  #R version 4.2.1 (2022-06-23 ucrt)
  #Platform: x86_64-w64-mingw32/x64 (64-bit)
  #Running under: Windows 10 x64 (build 19044)
  
  #Matrix products: default
  
  #locale:
  #[1] LC_COLLATE=Finnish_Finland.utf8  LC_CTYPE=Finnish_Finland.utf8    LC_MONETARY=Finnish_Finland.utf8
  #[4] LC_NUMERIC=C                     LC_TIME=Finnish_Finland.utf8    
  
  #attached base packages:
  #[1] stats     graphics  grDevices utils     datasets  methods   base     
  
  #other attached packages:
  #[1] lmeresampler_0.2.1.99999 lme4_1.1-29              Matrix_1.4-1            
  
  #loaded via a namespace (and not attached):
  #[1] Rcpp_1.0.8.3         lubridate_1.8.0      lattice_0.20-45      tidyr_1.2.1          prettyunits_1.1.1   
  #[6] ps_1.7.1             rprojroot_2.0.3      digest_0.6.29        utf8_1.2.2           R6_2.5.1            
  #[11] plyr_1.8.7           HLMdiag_0.5.0        ggplot2_3.3.6        pillar_1.8.1         rlang_1.0.6         
  #[16] curl_4.3.2           rstudioapi_0.13      minqa_1.2.4          callr_3.7.0          nloptr_2.0.3        
  #[21] desc_1.4.1           devtools_2.4.3       diagonals_6.4.0      splines_4.2.1        statmod_1.4.37      
  #[26] stringr_1.4.1        munsell_0.5.0        compiler_4.2.1       janitor_2.1.0        pkgconfig_2.0.3     
  #[31] pkgbuild_1.3.1       mgcv_1.8-40          tidyselect_1.1.2     tibble_3.1.7         fansi_1.0.3         
  #[36] crayon_1.5.1         dplyr_1.0.10         withr_2.5.0          MASS_7.3-57          distributional_0.3.1
  #[41] ggdist_3.2.0         grid_4.2.1           nlme_3.1-158         gtable_0.3.1         lifecycle_1.0.3     
  #[46] magrittr_2.0.3       scales_1.2.1         cli_3.3.0            stringi_1.7.8        cachem_1.0.6        
  #[51] farver_2.1.1         reshape2_1.4.4       fs_1.5.2             remotes_2.4.2        snakecase_0.11.0    
  #[56] ellipsis_0.3.2       generics_0.1.3       vctrs_0.4.1          boot_1.3-28          tools_4.2.1         
  #[61] forcats_0.5.2        rcmdcheck_1.4.0      glue_1.6.2           purrr_0.3.4          processx_3.6.1      
  #[66] pkgload_1.3.0        fastmap_1.1.0        nlmeU_0.70-9         colorspace_2.0-3     xopen_1.0.0         
  #[71] sessioninfo_1.2.2    memoise_2.0.1        usethis_2.1.6          
  
  #RStudio 2022.02.0+443 "Prairie Trillium" Release (9f7969398b90468440a501cf065295d9050bb776, 2022-02-16) for Windows
  #Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.8 Chrome/69.0.3497.128 Safari/537.36
  
  #Obtained with the following code:
  
  #library(lme4)
  #library(lmeresampler)
  #model <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
  #nsim <- 100
  #set.seed(123)
  #B071022ref <- bootstrap(model, .f = fixef, type = "residual", B = nsim, rbootnoise = 0.0001)
  
  #About Ubuntu R-CMD-checks. It was noticed that the results produced with Linux systems have insignificant
  #deviations from the results of other systems. Small relative deviations around 10^-10 were observed between
  #the reference data created using Windows and test data created using Linux. Thus, exactly the same
  #reproduction of the reference data is not expected. Maximum relative deviations of 0.001 between the
  #results created using different operating systems is allowed. The linux-test data (not used in this test_that())
  #on which the above is based on was obtained with the following specs:

  #> sessionInfo()
  #R version 4.2.1 (2022-06-23)
  #Platform: x86_64-pc-linux-gnu (64-bit)
  #Running under: CentOS Linux 7 (Core)

  #Matrix products: default
  #BLAS/LAPACK: /home/opt/easybuild/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib64/libflexiblas.so.3.2

  #locale:
  #[1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
  #[3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
  #[5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
  #[7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
  #[9] LC_ADDRESS=C               LC_TELEPHONE=C
  #[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

  #attached base packages:
  #[1] stats     graphics  grDevices utils     datasets  methods   base

  #other attached packages:
  #[1] lmeresampler_0.2.4 lme4_1.1-29        Matrix_1.4-1

  #loaded via a namespace (and not attached):
  #[1] Rcpp_1.0.8.3         plyr_1.8.7           pillar_1.7.0
  #[4] compiler_4.2.1       nloptr_2.0.3         forcats_0.5.1
  #[7] tools_4.2.1          statmod_1.4.36       boot_1.3-28
  #[10] lubridate_1.8.0      lifecycle_1.0.3      tibble_3.1.7
  #[13] nlme_3.1-158         gtable_0.3.0         lattice_0.20-45
  #[16] mgcv_1.8-40          pkgconfig_2.0.3      HLMdiag_0.5.0
  #[19] rlang_1.0.6          cli_3.6.0            DBI_1.1.3
  #[22] stringr_1.4.0        dplyr_1.0.9          janitor_2.2.0
  #[25] generics_0.1.2       vctrs_0.5.2          diagonals_6.4.0
  #[28] grid_4.2.1           tidyselect_1.1.2     nlmeU_0.70-9
  #[31] snakecase_0.11.0     glue_1.6.2           R6_2.5.1
  #[34] fansi_1.0.3          distributional_0.3.0 minqa_1.2.4
  #[37] tidyr_1.2.0          farver_2.1.0         reshape2_1.4.4
  #[40] ggplot2_3.4.1        purrr_0.3.4          magrittr_2.0.3
  #[43] scales_1.2.0         ellipsis_0.3.2       MASS_7.3-57
  #[46] splines_4.2.1        ggdist_3.2.1         assertthat_0.2.1
  #[49] colorspace_2.0-3     utf8_1.2.2           stringi_1.7.6
  #[52] munsell_0.5.0        crayon_1.5.1


  model <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
  nsim <- 100
  set.seed(123)
  boo <- bootstrap(model, .f = fixef, type = "residual", B = nsim, rbootnoise = 0.0001)
  
  maxreldev <- (B071022ref[["stats"]][,2:5] - boo[["stats"]][,2:5])/B071022ref[["stats"]][,2:5]
  maxreldev <- max(abs(maxreldev))
  comparison <- (maxreldev < 1)
  expect_true(comparison, info = NULL, label = NULL)

  maxreldev <- (B071022ref[["replicates"]] - boo[["replicates"]])/B071022ref[["replicates"]]
  maxreldev <- max(abs(maxreldev))
  comparison <- (maxreldev < 1)
  expect_true(comparison, info = NULL, label = NULL)

})

test_that("Capture possible absurdly large effects of rbootnoise = 0.0001 on rep.mean and se, an implication of major technical issues. Note, the same seed cannot be followed exactly due to the random noise the rbootnoise feature generates!",{

  model <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
  nsim <- 2000
  set.seed(123)
  booref <- bootstrap(model, .f = fixef, type = "residual", B = nsim)
  set.seed(123) #Note, the set.seed(123) will not be followed as above due to the additional random noise generation!
  boo <- bootstrap(model, .f = fixef, type = "residual", B = nsim, rbootnoise = 0.0001)
  
  expect_false(all.equal(boo, booref))
  
  boodif <- (boo[["stats"]][["rep.mean"]] - booref[["stats"]][["rep.mean"]])/booref[["stats"]][["rep.mean"]]
  expect_true(max(abs(boodif)) < 100)
  
  boodif <- (boo[["stats"]][["se"]] - booref[["stats"]][["se"]])/booref[["stats"]][["se"]]
  expect_true((max(abs(boodif)) < 100))
  
})
