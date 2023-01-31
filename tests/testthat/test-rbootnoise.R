library(lme4, quietly = TRUE)

data(jsp728, package = "lmeresampler")

# ==============================================================================
context("rbootnoise")
# ==============================================================================

test_that("compare rbootnoise = 0 to lmeresampler 0.2.2 results before the implementation of the feature",{
  skip_on_cran()
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
  
  model <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
  nsim <- 100
  set.seed(123)
  boo <- bootstrap(model, .f = fixef, type = "residual", B = nsim)
  
  comparison <- identical(boo[["stats"]], A071022ref[["stats"]])
  expect_true(comparison, info = NULL, label = NULL)
  comparison <- identical(boo[["replicates"]], A071022ref[["replicates"]])
  expect_true(comparison, info = NULL, label = NULL)
})

test_that("compare rbootnoise = 0.0001 to lmeresampler 0.2.2 results before the implementation of the feature",{
  skip_on_cran()
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
  
  model <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
  nsim <- 100
  set.seed(123)
  boo <- bootstrap(model, .f = fixef, type = "residual", B = nsim, rbootnoise = 0.0001)
  
  comparison <- identical(boo[["stats"]], A071022ref[["stats"]])
  expect_false(comparison, info = NULL, label = NULL)
  comparison <- identical(boo[["replicates"]], A071022ref[["replicates"]])
  expect_false(comparison, info = NULL, label = NULL)
})

test_that("compare rbootnoise = 0.0001 to the results of the first implementation of the feature",{
  skip_on_cran()
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
  
  model <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
  nsim <- 100
  set.seed(123)
  boo <- bootstrap(model, .f = fixef, type = "residual", B = nsim, rbootnoise = 0.0001)
  
  comparison <- identical(boo[["stats"]], B071022ref[["stats"]])
  expect_true(comparison, info = NULL, label = NULL)
  comparison <- identical(boo[["replicates"]], B071022ref[["replicates"]])
  expect_true(comparison, info = NULL, label = NULL)
})

test_that("verify the small effect of rbootnoise = 0.0001 on rep.mean (<5%) and se (<1%)",{
  skip_on_cran()
  model <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
  nsim <- 2000
  set.seed(123)
  booref <- bootstrap(model, .f = fixef, type = "residual", B = nsim)
  set.seed(123) #Note, the set.seed(123) will not be followed as above due to the additional random noise generation!
  boo <- bootstrap(model, .f = fixef, type = "residual", B = nsim, rbootnoise = 0.0001)
  
  expect_false(identical(boo, booref))
  
  boodif <- (boo[["stats"]][["rep.mean"]] - booref[["stats"]][["rep.mean"]])/booref[["stats"]][["rep.mean"]]*100
  expect_true(max(abs(boodif)) < 5)
  
  boodif <- (boo[["stats"]][["se"]] - booref[["stats"]][["se"]])/booref[["stats"]][["se"]]*100
  expect_true((max(abs(boodif)) < 1))
  
})
