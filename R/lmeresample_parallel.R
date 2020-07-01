library(lmeresampler)
library(lme4)

# fit model
vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)

# run sequential parametric bootstrap
library(tictoc)
tic()
b_nopar  <- bootstrap(vcmodA, fn = fixef, type = "parametric", B = 5000)
toc()

# run sequential cases bootstrap
tic()
boo2 <- bootstrap(model = vcmodA, fn = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
toc()



# run sequential cgr bootstrap
tic()
boo4 <- bootstrap(model = vcmodA, fn = fixef, type = "cgr", B = 1000)
toc()

# run sequential resid bootstrap
tic()
boo4 <- bootstrap(model = vcmodA, fn = fixef, type = "residual", B = 1000)
toc()

# run sequential reb bootstrap
tic()
boo5 <- bootstrap(model = vcmodA, fn = fixef, type = "reb", B = 1000, reb_typ = 0)
toc()


# Using foreach to parallelize
library(foreach)
library(purrr)

combine <- function(...) {
  boot_list <- list(...)
  combo_stat <- map_dfr(boot_list, ~as.data.frame(.x$t))
  combo_r <- sum(map_dbl(boot_list, ~.x$R))
  RES <- boot_list[[1]]
  RES$t <- combo_stat
  RES$R <- combo_r
  return(RES)
}

library(doParallel)
registerDoParallel(cores = 2)

# parallel parametric bootstrap, 2 workers
tic()
b_parallel2 <- foreach(B = rep(2500, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(vcmodA, fn = fixef, type = "parametric", B = B)
toc()

# parallel cases bootstrap, 2 workers
# Why is the parallel version slower??!!
tic()
boo2_parallel <- foreach(B = rep(50, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(model = vcmodA, fn = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
toc()

profvis::profvis({
  boo2_parallel <- foreach(B = rep(50, 2), .combine = combine, .packages = "lmeresampler") %dopar%
    bootstrap(model = vcmodA, fn = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
})

# Let's investigate .cases.resamp() first...

repeated_case_resamp <- function(data, cluster, resample, B){
  purrr::map(1:B, ~.cases.resamp(dat = data, cluster = cluster, resample = resample))
}

# Change dopar to do for sequential
# For 2000 reps, sequential ~ 10.726s; parallel ~7.07 s (so it's about 1.5x faster)
tic()
resampled_cases <- foreach(B = rep(1000, 2), .combine = append, .packages = "lmeresampler") %dopar%
  repeated_case_resamp(data = vcmodA@frame, cluster = c(rev(names(lme4::getME(vcmodA, "flist"))), ".id"), resample = c(TRUE, FALSE), B = B)
toc()


# Now, let's see what happens if we add the model refit to the resampling bit...

.cases.resamp.refit <- function(model, dat, cluster, resample, fn) {
  form <- model@call$formula
  reml <- lme4::isREML(model)
  
  # exit early for trivial data
  if(nrow(dat) == 1 || all(resample==FALSE))
    return(dat)
  
  # ver <- as.numeric_version(packageVersion("dplyr"))
  res <- dat
  
  for(i in 1:length(cluster)) {
    
    if(i==1 & resample[i]) {
      dots <- cluster[1]
      grouped <- dplyr::group_by_(res, dots)
      g_rows <- dplyr::group_rows(grouped)
      # g_rows <- ifelse(ver >= "0.8.0", dplyr::group_rows(grouped), attributes(grouped)$indices)
      cls <- sample(seq_along(g_rows), replace = resample[i])
      idx <- unlist(g_rows[cls], recursive = FALSE)
      res <- res[idx, ]
    } else{
      if(i == length(cluster) & resample[i]) {
        dots <- cluster[-i]
        grouped <- dplyr::group_by_(res, .dots = dots)
        res <- dplyr::sample_frac(grouped, size = 1, replace = TRUE)
      } else{
        if(resample[i]) {
          dots <- cluster[i]
          res <- split(res, res[, cluster[1:(i-1)]], drop = TRUE)
          res <- plyr::ldply(res, function(df) {
            grouped <- dplyr::group_by_(df, .dots = dots)
            g_rows <- dplyr::group_rows(grouped)
            # g_rows <- ifelse(ver >= "0.8.0", dplyr::group_rows(grouped), attributes(grouped)$indices)
            cls <- sample(seq_along(g_rows), replace = resample[i])
            idx <- unlist(g_rows[cls], recursive = FALSE)
            grouped[idx, ]
          }, .id = NULL)
        }
      }
    }
    
    
  }
  # fn(update(model, newdata = res, REML = reml)) 
  fn(lme4::lmer(formula = form, data = res, REML = reml))
}

repeated_case_resamp_refit <- function(model, data, cluster, resample, fn, B){
  purrr::map(1:B, ~.cases.resamp.refit(model = model, dat = data, cluster = cluster, resample = resample, fn = fixef))
}

## Timings using lmer() for full refit
# For 2000 reps, sequential ~ 55.574s; parallel ~30.803 s (so it's about 1.8x faster)
# improvement was only about 1.5x for smaller B
## Timings using update() to refit model
# For 2000 reps, sequential ~ 62.258s; parallel ~33.231 s (so it's about 1.87x faster)
## Not sure whether to use update() or lmer()
tic()
bootstrap_stats <- foreach(B = rep(1000, 2), .combine = dplyr::bind_rows, .packages = "lmeresampler") %dopar%
  repeated_case_resamp_refit(
    model = vcmodA, 
    dat = vcmodA@frame, 
    cluster = c(rev(names(lme4::getME(vcmodA, "flist"))), ".id"), 
    resample = c(TRUE, FALSE), 
    fn = fixef,
    B = B
  )
toc()


# parallel residual bootstrap, 2 workers
tic()
boo4_parallel <- foreach(B = rep(500, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(model = vcmodA, fn = fixef, type = "residual", B = B)
toc()

# parallel cgr, 2 workers
tic()
boo3_parallel <- foreach(B = rep(500, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(model = vcmodA, fn = fixef, type = "cgr", B = B)
toc()

# parallel reb, 2 workers
tic()
boo5_parallel <- foreach(B = rep(500, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(model = vcmodA, fn = fixef, type = "reb", B = B, reb_typ = 0)
toc()

tic()
boo5_parallel <- foreach(B = rep(500, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(model = vcmodA, fn = fixef, type = "reb", B = B, reb_typ = 1)
toc()

