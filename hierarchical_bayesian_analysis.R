For the very curious: We can also model varying effects of parenting between educational groups. For this, we use Bayesian estimation in the `rstanarm`package, because estimation of hierarchical models is less reliable with other methods.

```{r}
library(rstanarm)
m = stan_lmer(ADHD ~
                poly(positive,2) + 
                poly(inconsistent,2) + 
                gender + birthmonth + 
                ( poly(positive,2) + 
                    poly(inconsistent,2) | mEDU),
              my_data, 
              QR= T,
              algorithm = "meanfield")

effects_i = data.frame(ggpredict(m,terms = c("inconsistent","mEDU"), type = "re"))
effects_p = data.frame(ggpredict(m,terms = c("positive","mEDU"), type = "re"))

plot_model(m, type = "re")

ggplot(effects_i, aes(x= x, y = predicted)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, x=x), alpha = 0.3) + 
  facet_wrap(~group)

ggplot(effects_p, aes(x= x, y = predicted)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, x=x), alpha = 0.3) + 
  facet_wrap(~group)

```

(Note we used only approximate Bayesian estimation here, because it would take to long with full Bayes).

With the package `brms` we could do this hierarchical regression _and_ use a beta binomial model. We keep this for another time.
