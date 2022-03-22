#' @title Bootstrap Nested Linear Mixed-Effects Models
#'
#' @description
#' Perform various bootstrap process for nested linear mixed effects (LMEs) models including:
#' parametric, residual, cases, wild, and REB bootstraps.
#'
#' @export
#' @param model The model object you wish to bootstrap.
#' @param .f A function returning the statistic(s) of interest.
#' @param type A character string indicating the type of bootstrap that is being
#'    requested. Possible values are \code{"parametric"}, \code{"residual"}, 
#'    \code{"case"}, \code{"wild"}, \code{"wild"}, or \code{"reb"} 
#'    (random effect block bootstrap).
#' @param B The number of bootstrap resamples.
#' @param resample A logical vector specifying whether each level of the model 
#'    should be resampled in the cases bootstrap. The levels should be specified 
#'    from the highest level (largest cluster) of the hierarchy to the lowest 
#'    (observation-level); for example for students within a school, specify the 
#'    school level first, then the student level.
#' @param reb_type Specification of what random effect block bootstrap version to
#'     implement. Possible values are \code{0}, \code{1} or \code{2}.
#' @param hccme either \code{"hc2"} or \code{"hc3"}, indicating which 
#'     heteroscedasticity consistent covariance matrix estimator to use.
#' @param aux.dist either \code{"f1"} or \code{"f2"} indicating which auxiliary 
#'     distribution to draw the errors from
#' 
#' @details
#' All of the below methods have been implemented for nested linear mixed-effects
#' models fit by \code{lmer} (i.e., an \code{lmerMod} object) and \code{lme} 
#' (i.e., an \code{lmerMod} object). Details of the bootstrap procedures can be found
#' in the help file for that specific function.
#'   
#' @return 
#' The returned value is an object of class "lmeresamp". This is a list
#' with the following elements:
#' 
#' \itemize{
#' \item \code{observed}: the estimated values for the model parameters
#' \item \code{model}: the fitted model object
#' \item \code{.f}: the function call
#' \item \code{replicates}: a \eqn{B \times p} data frame of bootstrap values for each of the p model parameters,
#' \item \code{stats}: a tibble containing the \]code{observed}, \code{rep.mean} (bootstrap mean), 
#' \code{se} (bootstrap standard error), and \code{bias} values for each model parameter,
#' \item \code{B}: the number of bootstrap resamples performed
#' \item \code{data}: the data with which the model was fit
#' \item \code{seed}: a vector of randomly generated seeds that are used by the bootstrap
#' \item \code{type}: the type of bootstrap executed
#' \item \code{call}: the call to \code{bootstrap()} that the user
#' \item \code{message}: a list of length \code{B} giving any messages generated during refitting. An entry will be \code{NULL} if no message was generated.
#' \item \code{warning}: a list of length \code{B} giving any warnings generated during refitting. An entry will be \code{NULL} if no message was generated.
#' \item \code{error}: a list of length \code{B} giving any errors generated during refitting. An entry will be \code{NULL} if no message was generated.
#' }
#' 
#' @seealso 
#' \itemize{
#'   \item \code{\link{parametric_bootstrap}}, \code{\link{resid_bootstrap}},
#'      \code{\link{case_bootstrap}}, \code{\link{reb_bootstrap}}, 
#'      \code{\link{wild_bootstrap}} for more details on a specific bootstrap.
#'   \item \code{\link[lme4]{bootMer}} in the \pkg{lme4} package for an 
#'      implementation of (semi-)parametric bootstrap for mixed models.
#' }
#' 
#' @examples 
#' library(lme4) 
#' vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
#' 
#' ## you can write your own function to return stats, or use something like 'fixef'
#' mySumm <- function(.) { 
#'   s <- getME(., "sigma")
#'     c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
#' }
#'
#' ## running a parametric bootstrap 
#' set.seed(1234)
#' boo1 <- bootstrap(model = vcmodA, .f = mySumm, type = "parametric", B = 20)
#' 
#' \dontrun{
#' ## running a cases bootstrap - only resampling the schools
#' boo2 <- bootstrap(model = vcmodA, .f = mySumm, type = "case", B = 100, resample = c(TRUE, FALSE))
#' 
#' ## running a cases bootstrap - resampling the schools and students within the school
#' boo2 <- bootstrap(model = vcmodA, .f = mySumm, type = "case", B = 100, resample = c(TRUE, FALSE))
#' 
#' ## running a residual bootstrap
#' boo4 <- bootstrap(model = vcmodA, .f = mySumm, type = "residual", B = 100)
#' 
#' ## running an REB0 bootstrap
#' boo5 <- bootstrap(model = vcmodA, .f = mySumm, type = "reb", B = 100, reb_typ = 0)
#' }
#' 
#' ## to print results in a formatted way
#' print(boo1)
#' 
#'
#' @references
#'    Carpenter, J. R., Goldstein, H. and Rasbash, J. (2003) A novel bootstrap 
#'    procedure for assessing the relationship between class size and achievement. 
#'    \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, 
#'    \bold{52}, 431--443.
#'    
#'    Chambers, R. and Chandra, H. (2013) A random effect block bootstrap for 
#'    clustered data. \emph{Journal of Computational and Graphical Statistics}, 
#'    \bold{22}, 452--470.
#'    
#'    Morris, J. S. (2002) The BLUPs are not "best" when it comes to bootstrapping. 
#'    \emph{Statistics and Probability Letters}, \bold{56}, 425--430.
#'    
#'    Van der Leeden, R., Meijer, E. and Busing F. M. (2008) Resampling multilevel 
#'    models. In J. de Leeuw and E. Meijer, editors, \emph{Handbook of 
#'    Multilevel Analysis}, pages 401--433. New York: Springer.
#'    
#'    Bates, D., Maechler, M., Bolker, W., Walker, S. (2015).
#'    Fitting Linear Mixed-Effects Models Using lme4. \emph{Journal of
#'    Statistical Software}, \bold{67}, 1--48. doi:10.18637/jss.v067.i01.
#'    
#'    Modugno, L., & Giannerini, S. (2015). The Wild Bootstrap for 
#'    Multilevel Models. \emph{Communications in Statistics -- Theory and Methods}, 
#'    \bold{44}(22), 4812--4825.
bootstrap <- function(model, .f, type, B, resample = NULL, reb_type = NULL, hccme = NULL, aux.dist = NULL) {
  if(!type %in% c("parametric", "residual", "case", "wild", "reb"))
    stop("'type' must be one of 'parametric', 'residual', 'case', 'wild', or 'reb'")
  if(!is.null(reb_type))
    if(!reb_type %in% 0:2) 
      stop("'reb_type' must be either 0, 1, or 2")
  UseMethod("bootstrap", model)
}


#' @title Parametric Bootstrap for Nested LMEs
#'
#' @description
#' Generate parametric bootstrap replicates of a statistic for a nested linear 
#' mixed-effects model.
#' 
#' @details
#' The parametric bootstrap simulates bootstrap samples from the estimated 
#' distribution functions. That is, error terms and random effects are simulated
#' from their estimated normal distributions and are combined into bootstrap
#' samples via the fitted model equation.
#'
#' @export
#' @inheritParams bootstrap
#' 
#' @return 
#' The returned value is an object of class "lmeresamp".
#' 
#' @seealso 
#' \itemize{
#'   \item Examples are given in \code{\link{bootstrap}}
#'   \item \code{\link{parametric_bootstrap}}, \code{\link{resid_bootstrap}},
#'      \code{\link{case_bootstrap}}, \code{\link{reb_bootstrap}}, 
#'      \code{\link{wild_bootstrap}} for more details on a specific bootstrap.
#'   \item \code{\link[lme4]{bootMer}} in the \pkg{lme4} package for an 
#'      implementation of (semi-)parametric bootstrap for mixed models.
#' }
#'
#' @references
#'    Chambers, R. and Chandra, H. (2013) A random effect block bootstrap for 
#'    clustered data. \emph{Journal of Computational and Graphical Statistics}, 
#'    \bold{22}, 452--470.
#'    
#'    Van der Leeden, R., Meijer, E. and Busing F. M. (2008) Resampling multilevel 
#'    models. In J. de Leeuw and E. Meijer, editors, \emph{Handbook of 
#'    Multilevel Analysis}, pages 401--433. New York: Springer.
parametric_bootstrap <- function(model, .f, B) {
  UseMethod("parametric_bootstrap", model)
}


#' @title Cases Bootstrap for Nested LMEs
#'
#' @description
#' Generate cases bootstrap replicates of a statistic for a nested linear 
#' mixed-effects model.
#'
#' @details 
#' The cases bootstrap is a fully nonparametric bootstrap that resamples the data
#' with respect to the clusters in order to generate bootstrap samples. Depending 
#' on the nature of the data, the resampling can be done only for the higher-level 
#' cluster(s), only at the observation-level within a cluster, or at all levels.
#' See Van der Leeden et al. (2008) for a nice discussion of this decision. 
#' 
#' To resample a given level of the model, the corresponding entry in the logical 
#' vector specified in the \code{resample} parameter must be set to true. A few
#' examples are given below in terms of a two-level model where students are
#' clustered within schools:
#'
#' \itemize{
#'   \item To resample only the schools, set \code{resample = c(TRUE, FALSE)}.
#'   \item To resample only the students, set \code{resample = c(FALSE, TRUE)}.
#'   \item To resample both the students and the schools, set \code{resample = c(TRUE, TRUE)}.
#' }
#' 
#' 
#' @export
#' @inheritParams bootstrap
#'
#' @return 
#' The returned value is an object of class "lmeresamp".
#' 
#' @seealso 
#' \itemize{
#'   \item Examples are given in \code{\link{bootstrap}}
#'   \item \code{\link{parametric_bootstrap}}, \code{\link{resid_bootstrap}},
#'      \code{\link{case_bootstrap}}, \code{\link{reb_bootstrap}}, 
#'      \code{\link{wild_bootstrap}} for more details on a specific bootstrap.
#'   \item \code{\link[lme4]{bootMer}} in the \pkg{lme4} package for an 
#'      implementation of (semi-)parametric bootstrap for mixed models.
#' }
#'
#' @references
#'    Van der Leeden, R., Meijer, E. and Busing F. M. (2008) Resampling multilevel 
#'    models. In J. de Leeuw and E. Meijer, editors, \emph{Handbook of 
#'    Multilevel Analysis}, pages 401--433. New York: Springer.
case_bootstrap <- function(model, .f, B, resample) {
  UseMethod("case_bootstrap", model)
}

#' Residual Bootstrap for Nested LMEs
#'
#' @description
#' Generate semi-parametric residual bootstrap replicates of a statistic for a nested 
#' linear mixed-effects model.
#'
#' @export
#' @inheritParams bootstrap
#' 
#' @details 
#' The semi-parametric bootstrap algorithm implemented was outlined by  Carpenter,  
#' Goldstein and Rasbash (2003), and is referred to as the CGR bootstrap by some. 
#' The algorithm is outlined below:
#' \enumerate{
#'   \item Obtain the parameter estimates from the fitted model and calculate
#'      the estimated error terms and EBLUPs.
#'   \item Center and rescale the error terms and EBLUPs so that the empirical variance of
#'      these quantities is equal to estimated variance components from the model.
#'   \item Sample independently with replacement from the rescaled estimated error 
#'      terms and rescaled EBLUPs.
#'   \item Obtain bootstrap samples by combining the samples via the fitted model equation.
#'   \item Refit the model and extract the statistic(s) of interest.
#'   \item Repeat steps 3-5 B times.
#' }
#'
#' @return 
#' The returned value is an object of class "lmeresamp".
#' 
#' @seealso 
#' \itemize{
#'   \item Examples are given in \code{\link{bootstrap}}
#'   \item \code{\link{parametric_bootstrap}}, \code{\link{resid_bootstrap}},
#'      \code{\link{case_bootstrap}}, \code{\link{reb_bootstrap}}, 
#'      \code{\link{wild_bootstrap}} for more details on a specific bootstrap.
#'   \item \code{\link[lme4]{bootMer}} in the \pkg{lme4} package for an 
#'      implementation of (semi-)parametric bootstrap for mixed models.
#' }
#'
#' @references
#'    Carpenter, J. R., Goldstein, H. and Rasbash, J. (2003) A novel bootstrap 
#'    procedure for assessing the relationship between class size and achievement. 
#'    \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, 
#'    \bold{52}, 431--443.
resid_bootstrap <- function(model, .f, B) {
  UseMethod("resid_bootstrap", model)
}

#' @title REB Bootstrap for Two-Level Nested LMEs
#'
#' @description
#' Generate random effect block (REB) bootstrap replicates of a statistic for a 
#' two-level nested linear mixed-effects model.
#'
#' @details
#' The random effects block (REB) bootstrap was outlined by Chambers and Chandra (2013)
#' and has been developed for two-level nested linear mixed-effects (LME) models. 
#' Consider a two-level LME of the form
#' \deqn{y = X \beta + Z b + \epsilon}
#' 
#' The REB bootstrap algorithm (\code{type = 0}) is as follows:
#' \enumerate{
#'   \item Calculate the nonparametric residual quantities for the fitted model
#'   \itemize{
#'      \item marginal residuals \eqn{r = y - X\beta}
#'      \item predicted random effects \eqn{\tilde{b} = (Z^\prime Z)^{-1} Z^\prime r}
#'      \item error terms \eqn{\tilde{e} = r - Z \tilde{b}}
#'   }
#'   \item Take a simple random sample, with replacement, of the predicted random effects, \eqn{\tilde{b}}.
#'   \item Draw a simple random sample, with replacement, of the group (cluster) IDs.
#'       For each sampled cluster, draw a random sample, with replacement, of size 
#'       \eqn{n_i} from that cluster's vector of error terms, \eqn{\tilde{e}}.
#'   \item Generate bootstrap samples via the fitted model equation 
#'      \eqn{y = X \widehat{\beta} + Z \tilde{b} + \tilde{e}}
#'   \item Refit the model and extract the statistic(s) of interest.
#'   \item Repeat steps 2-5 B times.
#' }
#' 
#' Variation 1 (\code{type = 1}): 
#'    The first variation of the REB bootstrap zero centers and rescales the 
#'    residual quantities prior to resampling.
#' 
#' Variation 2 (\code{type = 2}):
#'    The second variation of the REB bootstrap scales the estimates and centers
#'    the bootstrap distributions (i.e., adjusts for bias) after REB bootstrapping.
#' 
#' @export
#' @inheritParams bootstrap
#'
#' @return 
#' The returned value is an object of class "lmeresamp".
#' 
#' @seealso 
#' \itemize{
#'   \item Examples are given in \code{\link{bootstrap}}
#'   \item \code{\link{parametric_bootstrap}}, \code{\link{resid_bootstrap}},
#'      \code{\link{case_bootstrap}}, \code{\link{reb_bootstrap}}, 
#'      \code{\link{wild_bootstrap}} for more details on a specific bootstrap.
#'   \item \code{\link[lme4]{bootMer}} in the \pkg{lme4} package for an 
#'      implementation of (semi-)parametric bootstrap for mixed models.
#' }
#'
#' @references
#'    Chambers, R. and Chandra, H. (2013) A random effect block bootstrap for 
#'    clustered data. \emph{Journal of Computational and Graphical Statistics}, 
#'    \bold{22}, 452--470.
reb_bootstrap <- function(model, .f, B, reb_type) {
  UseMethod("reb_bootstrap", model)
}

#' Wild Bootstrap for LME models
#'
#' @description
#' Generate wild bootstrap replicates of a statistic for a 
#' linear mixed-effects model.
#'
#' @export
#' @inheritParams bootstrap
#' 
#' @details 
#' The wild bootstrap algorithm for LMEs implemented here was outlined by  
#' Modugno & Giannerini (2015). The algorithm is outlined below:
#' \enumerate{
#'   \item Draw a random sample equal to the number of groups (clusters) from
#'      an auxillary distribution with mean zero and unit variance. 
#'      Denote these as \eqn{w_1, \ldots, w_g}.
#'   \item Calculate the selected heteroscedasticity consistent matrix estimator for
#'      the marginal residuals, \eqn{\tilde{v}_i}
#'   \item Generate bootstrap responses using the fitted equation: 
#'      \eqn{y_i^* = X_i \beta + \tilde{v}_i w_j}
#'   \item Refit the model and extract the statistic(s) of interest.
#'   \item Repeat steps 2-4 B times.
#' }
#'
#' @return 
#' The returned value is an object of class "lmeresamp".
#' 
#' @seealso 
#' \itemize{
#'   \item Examples are given in \code{\link{bootstrap}}
#'   \item \code{\link{parametric_bootstrap}}, \code{\link{resid_bootstrap}},
#'      \code{\link{case_bootstrap}}, \code{\link{reb_bootstrap}}, 
#'      \code{\link{wild_bootstrap}} for more details on a specific bootstrap.
#'   \item \code{\link[lme4]{bootMer}} in the \pkg{lme4} package for an 
#'      implementation of (semi-)parametric bootstrap for mixed models.
#' }
#'
#' @references
#'    Modugno, L., & Giannerini, S. (2015). The Wild Bootstrap for 
#'    Multilevel Models. \emph{Communications in Statistics -- Theory and Methods}, 
#'    \bold{44}(22), 4812--4825.
wild_bootstrap <- function(model, .f, B, hccme, aux.dist) {
  UseMethod("wild_bootstrap", model)
}


#' Extract parameters from LME models
#'
#' @description
#' A utility function that extracts the fixed effects and variance component
#' estimates from a fitted \code{merMod} or \code{lme} object.
#' 
#' @return 
#' A named vector of parameters.
#' 
#' @export
#' @param model the fitted model object, either an \code{merMod} or \code{lme} object.
extract_parameters <- function(model) {
  UseMethod("extract_parameters", model)
}