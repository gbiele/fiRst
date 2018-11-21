library(margins)
par(mar=c(3,3,2,1), mgp=c(2,.5,0), tck=-.01)
pcf = function(m) {
  bci = cbind(coef(m),
              confint(m))
  y = 1:nrow(bci)
  par(mar = c(3,10,1,1))
  plot(bci[,1],y,
       xlim = range(bci),
       pch = 16,
       yaxt = "n",
       ylab = "",
       xlab = "coefficient value")
  segments(x0 = bci[,2], x1 = bci[,3], y0 = y)
  axis(2,at = y, labels = rownames(bci), las = 2)
  abline(v = 0, lty = 2, col = "red")
}

par(mfrow = c(2,2))
hist(residuals(ADHD_model))
pcf(ADHD_model)
par(mar=c(3,3,2,1), mgp=c(2,.5,0), tck=-.01)
cplot(ADHD_model_bb,
      "positive",
      what = "prediction",
      ylim = quantile(my_data$ADHD,c(.025,.975)))
cplot(ADHD_model,
      "inconsistent",
      what = "prediction",
      ylim = quantile(my_data$ADHD,c(.025,.975)))


To see how useful it can be to "program" in R, lets plot the results for all outcomes. First, we consolidate all plotting into one function:
  
  ```{r}
plot_regression_results = function(m,outcome) {
  par(mar=c(3,3,2,1), mgp=c(2,.5,0), tck=-.01)
  par(mfrow = c(2,2))
  hist(residuals(m))
  pcf(m)
  title(outcome)
  par(mar=c(3,3,2,1), mgp=c(2,.5,0), tck=-.01)
  cplot(m,
        "positive",
        what = "prediction",
        ylim = quantile(my_data[,outcome],c(.025,.975)))
  title(outcome)
  cplot(m,
        "inconsistent",
        what = "prediction",
        ylim = quantile(my_data[,outcome],c(.025,.975)))
  title(outcome)
}
```

```{r}
outcomes = c("OD","CD","SMFQ","CCC","SPRAAK")
for (o in outcomes) {
  reg_mod_str = paste(o,"~poly(positive,2) + poly(inconsistent,2)")
  lm_fit = lm(as.formula(reg_mod_str),
              data = my_data)
  plot_regression_results(lm_fit,o)
}
```

To quickly fo through all outcomes, we use the `paste` and `as.formula`functions to first glue together our regression formula as a string, and then convert it from string to a formula object.

```{r}
my_data = my_data[complete.cases(my_data[,c("mEDU","mAge","parity","gender")]),]
outcomes = c("OD","CD","SMFQ","CCC","SPRAAK")
my_data$parity_n = as.numeric(my_data$parity)
my_data$mEDU_n = as.numeric(my_data$mEDU)
for (o in outcomes) {
  reg_mod_str = paste(o,"~poly(positive,2) + poly(inconsistent,2) + poly(mEDU_n,2) + poly(mAge,2) + poly(parity_n,2) + gender")
  lm_fit = lm(as.formula(reg_mod_str),
              data = my_data)
  plot_regression_results(lm_fit,o)
}
```
