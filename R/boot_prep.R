# bootstrap CIs:

## normal t

## these are not in the right order

con <- data.frame(lme4::confint.merMod(model))
colnames(con)[1] <- "norm.t.lower"
norm.t.lower <- con[1]
colnames(con)[2] <- "norm.t.upper"
norm.t.upper <- con[2]

## boot t

### This will need to be changed for nlme because accessing coefficients is different

# this gets all of the means (estimates) for boot_t calculation
model.fits <- lme4::getME(model, "beta")
model.fits <- append(model.fits, lme4::getME(model, "sigma"))
model.fits <- append(model.fits, lme4::getME(model, "theta"))

# sd of of estimates for boot_t calculation
out <- summary(model)
model.sds <- out$coefficients[,2] # fixef
model.sds <- append(model.sds, 0) # not sure what sd for sigma should be
model.sds <- append(model.sds, merTools::REsdExtract(model)) # random effects

# table of estimates and sds for boot_t calculation
t.stats <- cbind(model.fits, model.sds)
row.names(t.stats) <- colnames(replicates)

## andy's formula
boot_t <- (mean(x) - model@beta[x])*(sd(x)/sqrt(B))
model@beta - quantile(boot_t,c(0.975, 0.025)) * sd(model@beta)/sqrt(B)

# fix this once the math works out
boot.t.lower <- mapply(replicates[1:nrow(tstar)], model.fits, function(x) {
  model.fit[x]-qt(0.975,df=B-1)*sd(x)/sqrt(B-1), 8))
})

boot.t.upper <- apply(replicates[1:nrow(tstar)], 2, function(x) {
  round(mean(x)-qt(0.025,df = B-1)*sd(x)/sqrt(B), 8)
})

## percentile t

perc.t.lower <- apply(replicates, 2, function(x) {
  round(quantile(x, 0.025), 8)
})

perc.t.upper <- apply(replicates, 2, function(x) {
  round(quantile(x, 0.975), 8)
})

# put them all together 
norm.t.cis <-  data.frame(cbind(norm.t.lower, norm.t.upper))
other.cis <- data.frame(cbind(boot.t.lower, boot.t.upper, perc.t.lower, perc.t.upper))
cat(paste("Normal 95% t-interval 95%"))
print(norm.t.cis)
cat(paste("95% bootstrap-t and percentile confidence intervals:"))
print(other.cis)
