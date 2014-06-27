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