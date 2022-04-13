#' @title Calculate Bootstrap p-values for fixed effects
#' 
bootstrap_pvals.merMod <- function(model, type, B, resample = NULL,
                                   reb_type = NULL, 
                                   hccme = NULL, 
                                   aux.dist = NULL) {
  
  mod_summary <- summary(model)
  tstats <- stats::coef(mod_summary)[, "t value"]
  trms <- names(tstats)
  xmat <- lme4::getME(model, "X")
  
  
  pvals <- rep(NA, length(tstats))
  for(i in seq_along(pvals)) {
    trm <- trms[i]
    if(trm == "(Intercept)") null_mod <- update(model, . ~ . - 1)
    else null_mod <- update(model, as.formula(paste(". ~ . -", trm)))
    
    if(type == "case") {
      gen_data <- bootstrap(null_mod, type = type, B = B, .refit = FALSE, 
                            resample = resample, orig_data = model@frame)
    } else {
      gen_data <- bootstrap(null_mod, type = type, B = B, .refit = FALSE, resample = resample,
                            reb_type = reb_type, hccme = hccme, aux.dist = aux.dist)
    }
    
    extract_t <- function(x) stats::coef(summary(x))[i, "t value"]
    
    if(type == "case") {
      boot_ts <- purrr::map_dbl(gen_data, ~refit_to_newdf(model, newdata = .x, .f = extract_t))
    } else{
      boot_ts <- unlist(refit_merMod(gen_data, model, .f = extract_t)$tstar)
    }
    
    pvals[i] <- (sum(abs(boot_ts) >= abs(tstats[i])) + 1) / (B + 1)
  }
  
  cbind(stats::coef(mod_summary), p.value = pvals)
  
}
