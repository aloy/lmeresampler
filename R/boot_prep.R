# bootstrap CIs:

## normal t

con <- as.data.frame(confint(model, level = 0.95))
colnames(con)[1] <- "norm.t.lower"
norm.t.lower <- con[1]
colnames(con)[2] <- "norm.t.upper"
norm.t.upper <- con[2]

## boot t

boot.t.lower <- apply(replicates[1:nrow(tstar)], 2, function(x) {
  round(mean(x)-qt(0.975,df=B-1)*sd(x)/sqrt(B-1), 8)
})

boot.t.upper <- apply(replicates[1:nrow(tstar)], 2, function(x) {
  round(mean(x)-qt(0.025,df = B-1)*sd(x)/sqrt(B), 8)
})

## percentile t


# put them all together 

cis <- as.data.frame(cbind(norm.t.lower, norm.t.upper, boot.t.lower, boot.t.upper))

cat("95% confidence intervals:", cis)
