pkgname <- "varTestnlme"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('varTestnlme')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("varCompTest")
### * varCompTest

flush(stderr()); flush(stdout())

### Name: varCompTest
### Title: Variance component testing
### Aliases: varCompTest varCompTest.lme varCompTest.lme4
###   varCompTest.saemix varCompTest.merMod

### ** Examples

# load lme4 package and example dataset
library(lme4)
data(Orthodont, package = "nlme")

# fit the two models under H1 and H0
m1 <- lmer(distance ~ 1 + Sex + age + age*Sex + 
(0 + age | Subject), data = Orthodont, REML = FALSE)
m0 <- lm(distance ~ 1 + Sex + age + age*Sex, data = Orthodont)

# compare them (order is important: m1 comes first)
varCompTest(m1,m0,pval.comp="bounds")

# using nlme
library(nlme)
m1 <- lme(distance ~ 1 + Sex + age + age*Sex, 
random = pdSymm(Subject ~ 1 + age), data = Orthodont, method = "ML")
m0 <- lme(distance ~ 1 + Sex, random = ~ 1 | Subject, data = Orthodont, method = "ML")

varCompTest(m1,m0)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
