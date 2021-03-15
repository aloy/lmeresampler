### Provided in lme4/predict.R
#' @importFrom stats family fitted model.frame model.response rbinom 
#' rgamma rnbinom rnorm rpois weights
gaussian_simfun <- function(object, nsim, ftd=fitted(object),
                            wts=weights(object)) {
  
  if (any(wts != 1)) warning("ignoring prior weights")
  rnorm(nsim*length(ftd), ftd, sd=sigma(object))
}

binomial_simfun <- function(object, nsim, ftd=fitted(object),
                            wts=weights(object)) {
  n <- length(ftd)
  ntot <- n*nsim
  if (any(wts %% 1 != 0))
    stop("cannot simulate from non-integer prior.weights")
  ## Try to figure out if the original data were
  ## proportions, a factor or a two-column matrix
  if (!is.null(m <- model.frame(object))) {
    y <- model.response(m)
    if(is.factor(y)) {
      ## ignore weights
      yy <- factor(1+rbinom(ntot, size = 1, prob = ftd),
                   labels = levels(y))
      split(yy, rep(seq_len(nsim), each = n))
    } else if(is.matrix(y) && ncol(y) == 2) {
      yy <- vector("list", nsim)
      for (i in seq_len(nsim)) {
        Y <- rbinom(n, size = wts, prob = ftd)
        YY <- cbind(Y, wts - Y)
        colnames(YY) <- colnames(y)
        yy[[i]] <- YY
      }
      yy
    } else
      rbinom(ntot, size = wts, prob = ftd)/wts
  } else rbinom(ntot, size = wts, prob = ftd)/wts
}

poisson_simfun <- function(object, nsim, ftd=fitted(object),
                           wts=weights(object)) {
  ## A Poisson GLM has dispersion fixed at 1, so prior weights
  ## do not have a simple unambiguous interpretation:
  ## they might be frequency weights or indicate averages.
  wts <- weights(object)
  if (any(wts != 1)) warning("ignoring prior weights")
  rpois(nsim*length(ftd), ftd)
}


## FIXME: need a gamma.shape.merMod method in order for this to work.
##        (see initial shot at gamma.shape.merMod below)
Gamma_simfun <- function(object, nsim, ftd=fitted(object),
                         wts=weights(object)) {
  if (any(wts != 1)) message("using weights to scale shape parameter")
  ## used to use gamma.shape(), but sigma() is more general
  ## (wouldn't work *outside* of the merMod context though)
  shape <- sigma(object)*wts
  rgamma(nsim*length(ftd), shape = shape, rate = shape/ftd)
}

gamma.shape.merMod <- function(object, ...) {
  if(family(object)$family != "Gamma")
    stop("Can not fit gamma shape parameter because Gamma family not used")
  
  y <- lme4::getME(object, "y")
  mu <- lme4::getME(object, "mu")
  w <- weights(object)
  # Sec 8.3.2 (MN)
  L <- w*(log(y/mu)-((y-mu)/mu))
  dev <- -2*sum(L)
  # Eqs. between 8.2 & 8.3 (MN)
  Dbar <- dev/length(y)
  structure(list(alpha = (6+2*Dbar)/(Dbar*(6+Dbar)),
                 SE = NA), # FIXME: obtain standard error
            class = "gamma.shape")
}

#' @importFrom statmod rinvgauss
inverse.gaussian_simfun <- function(object, nsim, ftd=fitted(object),
                                    wts = weights(object)) {
  if (any(wts != 1)) message("using weights as inverse variances")
  statmod::rinvgauss(nsim * length(ftd), mean = ftd,
                     shape= wts/sigma(object))
}

## in the original MASS version, .Theta is assigned into the environment
## (triggers a NOTE in R CMD check)
## modified from @aosmith16 GH contribution

negative.binomial_simfun <- function (object, nsim,
                                      ftd = fitted(object),
                                      wts=weights(object))
{
  
  if (any(wts != 1))
    warning("ignoring prior weights")
  theta <- lme4::getME(object, "glmer.nb.theta")
  rnbinom(nsim * length(ftd), mu = ftd, size = theta)
}

simfunList <- list(gaussian = gaussian_simfun,
                   binomial = binomial_simfun,
                   poisson  = poisson_simfun,
                   Gamma    = Gamma_simfun,
                   negative.binomial = negative.binomial_simfun,
                   inverse.gaussian = inverse.gaussian_simfun)