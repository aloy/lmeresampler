library(lme4)
library(nlme)

boot <- function (model, fn, type){
  switch(type,
         par = parametric(model, fn),
         res = residual(model, fn),
         case = case(model, fn),
         cgr = cgr(model, fn),
         reb = reb(model, fn, reb_type = 0),
         reb1 = reb(model, fn, reb_type = 1),
         reb2 = reb(model, fn, reb_type = 2))
}

parametric <- function (model, fn){
  
}

residual <- function (model, fn){
  
}

case <- function (model, fn){
  
}

cgr <- function (model, fn){
  
}

reb <- function (model, fn, reb_type){
  if(reb_type = 1){
    # Call reb1 here
  }
  # reb code here
  if(reb_type = 2){
    # Call reb2 here
  }
}

reb1 <- function (model, fn){
  
}

reb2 <- function (model, fn){
  
}