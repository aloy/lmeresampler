#' @title Calculate Bootstrap p-values for fixed effects
#' 
#' @inheritParams bootstrap
#' @export
bootstrap_pvals <- function(model, type, B, resample = NULL,
                            reb_type = NULL, 
                            hccme = NULL, 
                            aux.dist = NULL) {
  if(!type %in% c("parametric", "residual", "case", "wild", "reb"))
    stop("'type' must be one of 'parametric', 'residual', 'case', 'wild', or 'reb'")
  if(!is.null(reb_type))
    if(!reb_type %in% 0:2) 
      stop("'reb_type' must be either 0, 1, or 2")
  UseMethod("bootstrap_pvals", model)

}  

#' @rdname bootstrap_pvals
#' @export
#' @method bootstrap_pvals merMod
bootstrap_pvals.merMod <- function(model, type, B, resample = NULL,
                                   reb_type = NULL, 
                                   hccme = NULL, 
                                   aux.dist = NULL) {
  
  mod_summary <- summary(model)
  tstats <- stats::coef(mod_summary)[, "t value"]
  trms <- names(tstats)
  # xmat <- lme4::getME(model, "X")
  
  
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
  
  coef_tbl <- cbind(stats::coef(mod_summary), p.value = pvals) %>%
    tibble::as_tibble(rownames = "term")
  
  structure(
    list(coefficients = coef_tbl, B = B, seed = .Random.seed, type = type),
    class = "coef_tbl"
  )
  
  
}

#' @rdname bootstrap_pvals
#' @export
#' @method bootstrap_pvals lme
bootstrap_pvals.lme <- function(model, type, B, resample = NULL,
                                   reb_type = NULL, 
                                   hccme = NULL, 
                                   aux.dist = NULL) {
  
  mod_summary <- summary(model)
  tstats <- stats::coef(mod_summary)[, "t-value"]
  trms <- names(tstats)
  
  
  pvals <- rep(NA, length(tstats))
  for(i in seq_along(pvals)) {
    trm <- trms[i]
    if(trm == "(Intercept)") null_mod <- update(model, . ~ . - 1)
    else null_mod <- update(model, as.formula(paste(". ~ . -", trm)))
    
    if(type == "case") {
      gen_data <- bootstrap(null_mod, type = type, B = B, .refit = FALSE, 
                            resample = resample)
    } else {
      gen_data <- bootstrap(null_mod, type = type, B = B, .refit = FALSE, resample = resample,
                            reb_type = reb_type, hccme = hccme, aux.dist = aux.dist)
    }
    
    extract_t <- function(x) stats::coef(summary(x))[i, "t-value"]
    
    if(type == "case") {
      boot_ts <- purrr::map(gen_data, ~updated.model(model, new.data = .x)) %>%
        purrr::map_dbl(~extract_t(.x))
    } else{
      boot_ts <- purrr::map(gen_data, ~updated.model(model, new.y = .x)) %>%
        purrr::map_dbl(~extract_t(.x))
    }
    
    pvals[i] <- (sum(abs(boot_ts) >= abs(tstats[i])) + 1) / (B + 1)
  }
  
  coef_tbl <- cbind(stats::coef(mod_summary), p.value = pvals) %>%
    tibble::as_tibble(rownames = "term")
  
  structure(
    list(coefficients = coef_tbl, B = B, seed = .Random.seed, type = type),
    class = "coef_tbl"
  )
  
}
