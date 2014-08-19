library(nlme)

data(sleepstudy, package = "lme4")

model <- lme(Reaction ~ Days, data = sleepstudy, random = ~Days|Subject)
fn <- fixef
B <- 10

### BEGIN PAR CODE ###
fn <- match.fun(fn)

model.fixef <- fixed.effects(model) # Extract fixed effects
for(i in 1:B){
  sim <- simulate.lme.data(model, nsim = 1, na.action = na.exclude)
  ystar[i] <- sim
}
ystar <- simulate.lme.data(model, nsim = 2, na.action = na.exclude)

t0 <- fn(model)


# Originally I had this: 
# tstar <- lapply(ystar, function(x) {
#   update(object = model, x ~ .)
#   fn(model)
# })
# But ran into this error:
# Error in eval(expr, envir, enclos) : object 'x' not found
# And when I change it to Reaction it runs but I do not think it works
# Refit the model and apply 'fn' to it using lapply
tstar <- lapply(ystar, function(x) {
  model.update <- update(object = model, fixed = x ~ .)
  t.res <- fn(model.update)
  return(t.res)
})

tstar <- do.call("cbind", tstar) # Can these be nested?
rownames(tstar) <- names(t0)

RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                      seed = .Random.seed, statistic = fn,
                      sim = "parametric", call = match.call()),
                 class = "boot")




####### Singular test ######
model <- lme(Reaction ~ Days, data = sleepstudy, random = ~Days|Subject)
fn <- fixed.effects
B <- 1

### BEGIN PAR CODE ###
fn <- match.fun(fn)

model.fixef <- fixed.effects(model) # Extract fixed effects
ystar <- simulate.lme.data(model, nsim = B, na.action = na.exclude)

t0 <- fn(model)

# Refit the model and apply 'fn' to it using lapply
model.update <- update(object = model, ystar ~ .)
anova(model.update, model) # This should not run if the update worked
tstar <- fn(model.update)

tstar <- do.call("cbind", tstar) # Can these be nested?
rownames(tstar) <- names(t0)

RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                      seed = .Random.seed, statistic = fn,
                      sim = "parametric", call = match.call()),
                 class = "boot")









################ Temp Sim FUN #######
simulate.lme.data <- function (object, nsim = 1, seed = as.integer(runif(1, 0, .Machine$integer.max)), 
          m2, method = c("REML", "ML"), niterEM = c(40, 200), useGen = FALSE, 
          ...) 
{
  if (inherits(nsim, "lm") || inherits(nsim, "lme")) 
    stop("order of arguments in 'simulate.lme' has changed to conform with generic in R-2.2.0", 
         domain = NA)
  getResults1 <- function(conLin, nIter, pdClass, REML, ssq, 
                          p, pp1) {
    unlist(.C(mixed_combined, as.double(conLin$Xy), as.integer(unlist(conLin$dims)), 
              double(ssq), as.integer(nIter), as.integer(pdClass), 
              as.integer(REML), logLik = double(1), R0 = double(pp1), 
              lRSS = double(1), info = integer(1))[c("info", "logLik")])
  }
  getResults2 <- function(conLin, reSt, REML, control) {
    lmeSt <- lmeStruct(reStruct = reStruct(reSt, REML = REML))
    attr(lmeSt, "conLin") <- conLin
    lmeSt <- Initialize(lmeSt, data = NULL, groups = NULL, 
                        control = control)
    attr(lmeSt, "conLin") <- MEdecomp(attr(lmeSt, "conLin"))
    aMs <- nlminb(c(coef(lmeSt)), function(lmePars) -logLik(lmeSt, 
                                                            lmePars), control = list(iter.max = control$msMaxIter, 
                                                                                     eval.max = control$msMaxEval, trace = control$msVerbose))
    c(info = aMs$flags[1], logLik = -aMs$value)
  }
  if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1)
  RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(seed)
  if (inherits(object, "lme")) {
    fit1 <- object
    object <- as.list(object$call[-1])
  }
  else {
    object <- as.list(match.call(lme, substitute(object))[-1])
    fit1 <- do.call("lme", object)
  }
  if (length(fit1$modelStruct) > 1) {
    stop("models with \"corStruct\" and/or \"varFunc\" objects not allowed")
  }
  reSt1 <- fit1$modelStruct$reStruct
  condL1 <- do.call("createConLin", object)
  pdClass1 <- unlist(lapply(reSt1, data.class))
  pdClass1 <- match(pdClass1, c("pdSymm", "pdDiag", "pdIdent", 
                                "pdCompSymm", "pdLogChol"), 0) - 1
  control1 <- lmeControl()
  if (!is.null(object$control)) {
    control1[names(object$control)] <- object$control
  }
  control1$niterEM <- niterEM[1]
  sig <- fit1$sigma
  DeltaInv <- pdMatrix(reSt1, factor = TRUE)
  for (i in names(DeltaInv)) {
    DeltaInv[[i]] <- sig * DeltaInv[[i]]
  }
  if (missing(useGen)) {
    useGen <- any(pdClass1 == -1)
  }
  nullD <- condL1$dims
  N <- nullD$N
  Q <- nullD$Q
  p1 <- nullD$ncol[Q + 1]
  pp11 <- p1 * (p1 + 1)
  ycol1 <- sum(nullD$ncol)
  qvec <- nullD$qvec[1:Q]
  ssq1 <- sum(qvec^2)
  csq1 <- cumsum(c(1, qvec[-Q]))
  csq2 <- cumsum(qvec)
  ngrp <- nullD$ngrps
  ind <- vector("list", Q)
  base <- condL1$Xy[, ycol1 - (nullD$ncol[Q + 1]:1), drop = FALSE] %*% 
    fixef(fit1)
  for (i in 1:Q) {
    ind[[i]] <- rep(1:ngrp[i], nullD$ZXlen[[i]])
  }
  value <- list(null = list())
  if (ML <- !is.na(match("ML", method))) {
    value$null$ML <- array(0, c(nsim, 2), list(1:nsim, c("info", 
                                                         "logLik")))
  }
  if (REML <- !is.na(match("REML", method))) {
    value$null$REML <- array(0, c(nsim, 2), list(1:nsim, 
                                                 c("info", "logLik")))
  }
  attr(value, "call") <- match.call()
  attr(value, "seed") <- seed
  ALT <- FALSE
  if (!missing(m2)) {
    ALT <- TRUE
    if (inherits(m2, "lme")) {
      fit2 <- m2
      m2 <- as.list(m2$call[-1])
    }
    else {
      m2 <- as.list(match.call(lme, substitute(m2))[-1])
      if (is.null(m2$random)) {
        m2$random <- asOneSidedFormula(object$fixed[-2])
      }
      aux <- object
      aux[names(m2)] <- m2
      m2 <- aux
      fit2 <- do.call("lme", m2)
    }
    if (length(fit2$modelStruct) > 1) {
      stop("models with \"corStruct\" and/or \"varFunc\" objects not allowed")
    }
    condL2 <- do.call("createConLin", m2)
    reSt2 <- fit2$modelStruct$reStruct
    control2 <- lmeControl()
    if (!is.null(m2$control)) {
      control2[names(m2$control)] <- m2$control
    }
    control2$niterEM <- niterEM[2]
    pdClass2 <- unlist(lapply(fit2$modelStruct$reStruct, 
                              data.class))
    pdClass2 <- match(pdClass2, c("pdSymm", "pdDiag", "pdIdent", 
                                  "pdCompSymm", "pdLogChol"), 0) - 1
    useGen <- useGen || any(pdClass2 == -1)
    altD <- condL2$dims
    ssq2 <- sum((altD$qvec[1:altD$Q])^2)
    p2 <- altD$ncol[altD$Q + 1]
    pp12 <- p2 * (p2 + 1)
    ycol2 <- sum(altD$ncol)
    if (ML) {
      value$alt$ML <- value$null$ML
    }
    if (REML) {
      value$alt$REML <- value$null$REML
    }
  }
  for (i in 1:nsim) {
    base2 <- base + rnorm(N, sd = sig)
    for (j in 1:Q) {
      base2 <- base2 + ((array(rnorm(ngrp[j] * qvec[j]), 
                               c(ngrp[j], qvec[j]), list(1:ngrp[j], NULL)) %*% 
                           DeltaInv[[j]])[ind[[j]], , drop = FALSE] * condL1$Xy[, 
                                                                                csq1[j]:csq2[j], drop = FALSE]) %*% rep(1, qvec[j])
    }
    condL1$Xy[, ycol1] <- base2
    if (REML) {
      if (useGen) {
        value$null$REML[i, ] <- getResults2(condL1, reSt1, 
                                            TRUE, control1)
      }
      else {
        value$null$REML[i, ] <- getResults1(condL1, niterEM[1], 
                                            pdClass1, TRUE, ssq1, p1, pp11)
      }
    }
    if (ML) {
      if (useGen) {
        value$null$ML[i, ] <- getResults2(condL1, reSt1, 
                                          FALSE, control1)
      }
      else {
        value$null$ML[i, ] <- getResults1(condL1, niterEM[1], 
                                          pdClass1, FALSE, ssq1, p1, pp11)
      }
    }
    if (ALT) {
      condL2$Xy[, ycol2] <- base2
      if (REML) {
        if (useGen) {
          value$alt$REML[i, ] <- getResults2(condL2, 
                                             reSt2, TRUE, control2)
        }
        else {
          value$alt$REML[i, ] <- getResults1(condL2, 
                                             niterEM[2], pdClass2, TRUE, ssq2, p2, pp12)
        }
      }
      if (ML) {
        if (useGen) {
          value$alt$ML[i, ] <- getResults2(condL2, reSt2, 
                                           FALSE, control2)
        }
        else {
          value$alt$ML[i, ] <- getResults1(condL2, niterEM[2], 
                                           pdClass2, FALSE, ssq2, p2, pp12)
        }
      }
    }
  }
  if (ML) {
    value$null$ML[, "logLik"] <- N * (log(N) - (1 + log(2 * 
                                                          pi)))/2 + value$null$ML[, "logLik"]
    if (ALT) {
      value$alt$ML[, "logLik"] <- N * (log(N) - (1 + log(2 * 
                                                           pi)))/2 + value$alt$ML[, "logLik"]
    }
  }
  if (REML) {
    value$null$REML[, "logLik"] <- (N - p1) * (log(N - p1) - 
                                                 (1 + log(2 * pi)))/2 + value$null$REML[, "logLik"]
    if (ALT) {
      value$alt$REML[, "logLik"] <- (N - p2) * (log(N - 
                                                      p2) - (1 + log(2 * pi)))/2 + value$alt$REML[, 
                                                                                                  "logLik"]
    }
  }
  attr(value, "df") <- p1 + length(coef(reSt1)) + 1
  if (ALT) {
    attr(value, "df") <- abs(attr(value, "df") - (p2 + length(coef(reSt2)) + 
                                                    1))
  }
  attr(value, "useGen") <- useGen
  class(value) <- "simulate.lme"
  assign(".Random.seed", RNGstate, envir = .GlobalEnv)
  value
  return(base2)
}
