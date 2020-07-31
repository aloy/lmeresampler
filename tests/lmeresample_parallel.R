library(lmeresampler)
library(lme4)
library(nlme)

# fit model
vcmodA <- lmer(mathAge11 ~ mathAge8 + mathAge8c + gender + class + (1 | school), data = jsp728)
vcmodB <- lme(mathAge11 ~ mathAge8 + gender + class, random = ~1|school, data = jsp728)

getME(vcmodA, "X")

## you can write your own function to return stats, or use something like 'fixef'
mySumm <- function(.) { 
  s <- lme4::getME(., "sigma")
  c(beta = lme4::getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
}

set.seed(1234)
# run sequential parametric bootstrap

b_resids  <- bootstrap(vcmodA, .f = residuals, type = "parametric", B = 100)

## lme4
library(tictoc)
tic()
b_nopar  <- bootstrap(vcmodA, .f = fixef, type = "parametric", B = 100)
toc()

## nlme
tic()
b_nopar2  <- bootstrap(vcmodB, .f = fixef, type = "parametric", B = 500)
toc()

# run sequential cases bootstrap

## lme4
tic()
boo2 <- bootstrap(model = vcmodA, .f = fixef, type = "case", B = 500, resample = c(TRUE, FALSE))
toc()

## nlme
tic()
boo2.2 <- bootstrap(model = vcmodB, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
toc()

# run sequential cgr bootstrap

## lme4
tic()
boo3 <- bootstrap(model = vcmodA, .f = fixef, type = "cgr", B = 100)
toc()

## nlme
tic()
boo3.2 <- bootstrap(model = vcmodB, .f = fixef, type = "cgr", B = 100)
toc()

# run sequential resid bootstrap
tic()
boo4 <- bootstrap(model = vcmodA, .f = fixef, type = "residual", B = 100, linked = TRUE)
toc()

## nlme
tic()
boo4.2 <- bootstrap(model = vcmodB, .f = fixef, type = "residual", B = 100, linked = TRUE)
toc()

# run sequential reb bootstrap
tic()
boo5 <- bootstrap(model = vcmodA, .f = fixef, type = "reb", B = 100, reb_type = 0)
toc()

tic()
boo5.2 <- bootstrap(model = vcmodB, .f = fixef, type = "reb", B = 100, reb_type = 0)
toc()


# Using foreach to parallelize
library(foreach)
library(purrr)

combine <- function(...) {
  boot_list <- list(...)
  combo_stat <- purrr::map_dfr(boot_list, ~as.data.frame(.x$t))
  combo_r <- sum(map_dbl(boot_list, ~.x$R))
  RES <- boot_list[[1]]
  RES$t <- combo_stat
  RES$R <- combo_r
  return(RES)
}

# ---------------------------------------------------------------------------------

# doParallel and foreach will automatically adapt to the user's system?

## CLUSTERING ##
# pro: system agnostic, con: slow
library(doParallel)
set.seed(1234)
cl <- snow::makeSOCKcluster(2)
doParallel::registerDoParallel(cl)

# parallel parametric bootstrap, 2 workers
# all packages being used need to explicitly specified in .packages
tic()
b_parallel2 <- foreach(B = rep(250, 2), .combine = combine, .packages = c("lmeresampler", "lme4")) %dopar%
  bootstrap(vcmodA, .f = fixef, type = "parametric", B = B)
toc()
stopCluster(cl)

# parallel cases bootstrap, 2 workers
# Why is the parallel version slower??!!
cl <- snow::makeSOCKcluster(2)
doParallel::registerDoParallel(cl)

tic()
boo2_parallel <- foreach(B = rep(250, 2), .combine = combine, .packages = c("lmeresampler", "lme4")) %dopar%
  bootstrap(model = vcmodA, .f = .fixef, type = "case", B = B, resample = c(TRUE, FALSE))
toc()
snow::stopCluster(cl)

profvis::profvis({
  boo2_parallel <- foreach(B = rep(50, 2), .combine = combine, .packages = "lmeresampler") %dopar%
    bootstrap(model = vcmodA, .f = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
})

snow::stopCluster(cl)

# ---------------------------------------------------------------------------------

## FORKING ##
# pro: faster, con: only works on UNIX systems (not Windows), so runtime will not be a Windows runtime

doParallel::registerDoParallel(cores = 2)

# parallel parametric bootstrap, 2 workers
tic()
b_parallel2 <- foreach(B = rep(250, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(vcmodA, .f = fixef, type = "parametric", B = B)
toc()

# parallel cases bootstrap, 2 workers
# Why is the parallel version slower??!!
tic()
boo2_parallel <- foreach(B = rep(250, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(model = vcmodA, .f = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
toc()

# ---------------------------------------------------------------------------------

# Let's investigate .cases.resamp() first...

repeated_case_resamp <- function(data, cluster, resample, B){
  purrr::map(1:B, ~.cases.resamp(dat = data, cluster = cluster, resample = resample))
}

snow::stopCluster(cl)

# Change dopar to do for sequential
# For 2000 reps, sequential ~ 10.726s; parallel ~7.07 s (so it's about 1.5x faster)
tic()
resampled_cases <- foreach(B = rep(1000, 2), .combine = append, .packages = "lmeresampler") %dopar%
  repeated_case_resamp(data = vcmodA@frame, cluster = c(rev(names(lme4::getME(vcmodA, "flist"))), ".id"), resample = c(TRUE, FALSE), B = B)
toc()

snow::stopCluster(cl)

# Now, let's see what happens if we add the model refit to the resampling bit...

.cases.resamp.refit <- function(model, dat, cluster, resample, .f) {
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
  # .f(update(model, newdata = res, REML = reml)) 
  .f(lme4::lmer(formula = form, data = res, REML = reml))
}

repeated_case_resamp_refit <- function(model, data, cluster, resample, .f, B){
  purrr::map(1:B, ~.cases.resamp.refit(model = model, dat = data, cluster = cluster, resample = resample, .f = fixef))
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
    .f = fixef,
    B = B
  )
toc()


# parallel residual bootstrap, 2 workers
tic()
boo4_parallel <- foreach(B = rep(500, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(model = vcmodA, .f = fixef, type = "residual", B = B)
toc()

# parallel cgr, 2 workers
tic()
boo3_parallel <- foreach(B = rep(500, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(model = vcmodA, .f = fixef, type = "cgr", B = B)
toc()

# parallel reb, 2 workers
tic()
boo5_parallel <- foreach(B = rep(500, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(model = vcmodA, .f = fixef, type = "reb", B = B, reb_typ = 0)
toc()

tic()
boo5_parallel <- foreach(B = rep(500, 2), .combine = combine, .packages = "lmeresampler") %dopar%
  bootstrap(model = vcmodA, .f = fixef, type = "reb", B = B, reb_typ = 1)
toc()
