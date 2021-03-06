library(lme4)
options(mc.cores = 4)

m = stan_lmer(ADHD ~
                poly(positive,2) + 
                poly(inconsistent,2) + 
                gender + birthmonth + 
                ( poly(positive,2) + 
                    poly(inconsistent,2) | mEDU),
              my_data, 
              algorithm = "meanfield",
              QR = T)

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
---
  title: "A first R analysis"
author: "Guido Biele"
date: "16 november 2018"
output:
  html_document: default
pdf_document: default
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE,
                      warning=FALSE,
                      message=FALSE)
```

## Introduction

This document aims to povide a shirt introduction into statistical analysis with R. The assumption is, that the reader is already familiar with basic concepts in R and the Rstuido IDE, and has some basic knowledge about statistics, for example from using SPSS. The goal here is to show how one can

* open an SPSS file
* merge seperate data sets
* visually inspect the data
* do a simple data reduction through an explatory factor analysis
* estimate the association between an exposure and an outcome with a linear regregssion
* visualize the results of this regression.

## Packages

The name R refers at the same time to 3 things:
  
  1. A programming language
2. A software for statistical analyses
3. An ecosystem of _packages_ that implement different statistical analyses.

R packages are written (typically by using R, the programming language) by researchers and available for free. These packages are crucial, because the capabilities of The basic R software are limited.

There are thousands of R packages. Therefore, CRAN [task view web sites](https://cran.rstudio.com/web/views/) are usful to give an overview about R packages for different application areas. CRAN stands for "Comprehensive R Archive Network".


For this introduction we will use following packes:
  
  * `haven` allows to read SPSS (or stata) data into R
* `VIM` for data imputation
* `psych` is a package for factor analysis and related things
* `ggplot` and `cowplot` for general plotting
* `sjPlot` and `ggeffect` for plotting of regression analysis results 
* `dplyr` for data wrangling
* `qualityTools` to make QQ plots for special models
* `apaTables` -- nomen est omen

Installing packages is easy!
  
  ```{r, eval=FALSE}
install.packages("haven")
install.packages("ggplot2")
install.packages("psych")
install.packages("sjplot")
install.packages("VIM")
```

One can also install multiple packages simulataneously as follows `install.packages(c("haven","psych"))`.

## Loading data

The problem we are looking at is the association between Parenting Style and child mental health symptoms. More specifically, we look at the association between self reported Parenting Style of well educated Norwegian mothers at child age 5 years and opposional behavior as well as mood and feeling at age 8. This data can be obtained from moba questionnaires 5 and 8.

If these data are available in SPSS (or stata) format, we can read them with the `haven` package, which we wirst have to load befor we can use it:
  
  ```{r}
library(haven)
library(VIM)
```

All packages comn with a documentation, and if one is lucky also with one or more _vignettes_, which are tutorial-style documents that show how the functions in the package can be used. This and more information about the packages can be found at https://cran.rstudio.com/web/packages/haven/. For other packages, just replace "haven" with the name of the package you are interested in.


OK, now we can load the data:
  ```{r, eval = F}
data_directory = "data/"
Q5aar = read_sav(paste0(data_directory,
                        "PDB439_Skjema5aar_v10.sav"))
Q5aar = data.frame(Q5aar)
```

So what have we done here?
  We first created a new variable `data_directory` and set it to the path to the directory where the data files are. Then, when loading the data into the data.frame `Q5aar`, we collated the path to the directory with the name of the data file by using the `paste0`command.

The reason for doing this is that the code becomes a bit easier to read (still not great ). For the same reason, I made a line break at the end of `Q5aar = read_sav(paste0(data_directory,`.

Note that if we would already be in the same directory as the data file, we could just have written `Q5aar = read_sav("PDB439_Skjema5aar_v10.sav")`. You can use the commands `getwd()` and `setwd()` to check and set you working directory (wd).

Next we select the variables we need. This are the pregancy ID and the child number, as well as the items for measuring parental style (we get SPSS variable names from [here](https://www.fhi.no/globalassets/dokumenterfiler/studier/moba/dokumenter/instrument-documentation-q-5year.pdf)). 

```{r, eval = F}
Q5aar = Q5aar[,c("PREG_ID_439","BARN_NR",
                 paste0("LL",526:534))]
names(Q5aar)
```

```{r, echo=F}
load("data/Q5aar.Rdata")
names(Q5aar)
```


Ok, these variable names are not very clear, lets rename them:
  
  ```{r}
names(Q5aar)[3:11] = paste0("PS.item.",1:9)
names(Q5aar)
```

It is important to understand the coding of the variables. The `haven` package makes SPSS' _Variable label_ and _Value label_ information available as attributes of the variables in the data.frame. We can access those as follows:

```{r}
attr(Q5aar$PS.item.1,"label")
attr(Q5aar$PS.item.1,"labels")
```

Here we see one anoying detail of SPSS Moba data: If the mother has ticked more than one box, the data are for our purposes missing, and not zero. We need to fix that.

To do this, we do a quick for loop:

```{r}
item_names = names(Q5aar)[-c(1,2)]
item_quest = c()
for (v in item_names) {
item_quest = c(item_quest,
strsplit(attr(Q5aar[,v],"label"),
split = ";")[[1]][2])
attr(Q5aar[, v],"labels") = attr(Q5aar[, v],"labels")[-1]
# need to use the "which", otherwise R stumbles over missing values
value_is_0 = which(Q5aar[,v] == 0)
Q5aar[value_is_0, v] = NA
}
table(as_factor(Q5aar$PS.item.1))
```

Missing data is always a problem, especially with questionnaire data from long lasting studies such as MoBa. Lets look at missing data by schowing a table that counts them:

```{r}
table_plus_missing = table(addNA(as_factor(Q5aar$PS.item.1)))
barplot(table_plus_missing)
table_plus_missing
```

There is lots of missing data, probably because MoBa changed the Parenting Style items throughout.

Lets look a little closer at the missing data:

```{r}
is_missing = is.na(Q5aar[,-c(1,2)])
hist(rowSums(is_missing))
```

Indeed there is a substantial number of cases for which als Parenting style itmes are missing. Let's remove those participants where more than 50% of parenting style items are missing, and impute the missing values for the remaining participants.  We use the function `hotdeck`from the package `VIM` to impute values
```{r}
Q5aar = Q5aar[rowMeans(is_missing[,item_names]) < .5,]
ps_imputations = hotdeck(as.matrix(Q5aar[,item_names]))
Q5aar[,item_names] = ps_imputations[,item_names]
```


## A quick exploratory factor analysis of parenting style items in MoBa

Our goal is to reduce the 9 items to a smaller set of variables describing parenting style. One way to do this is to use exploratory factor analysis (EFA) to learn the factor structure of the data. As a first step we'll check how many factors we need, using the `fa.parallel.poly` command from the `psych` package. 

As with many things in R, using up to date methods is not hard, as you'll see. The problem is more in knowing that these methods exist. If you do not know which method to use, try googling it. The challange is in finding out which of the many similar packages one should use. Here are  a few tipps:
  
  * google [R tutorial exploratory factor analysis](https://www.google.com/search?q=R+tutorial+exploratory+factor+analysis), 
* If a paper about a package is published in the [Journal of Statistical Software](https://www.jstatsoft.org), this is probably one of the better packages. dOn't worry about the name of the Journal! The articles typically include a tutorial about how to use the key functions.

Now lets move on with the exploratory factor analysis.
One thing we need to specify is that polychoric correlations are used, because we have ordinal data.

```{r}
library(psych)
describe(Q5aar[,item_names])
fa.parallel(Q5aar[,item_names], cor = "poly")
```

4 factors. I had hoped fewer, but lets continue and do an exploratory factor analysis.

```{r}
fa_ps = fa(Q5aar[,item_names],
nfactors = 4, 
cor = "poly")
fa_ps
```

We can at first look at the RMSEA. The value is 0.06, which is not great but OK.

This are the factor loadings: 

```{r}
plot(fa_ps)
fa.diagram(fa_ps)
```

This is not very clear, because a number of items (especially item 3) have high loadings for multiple factors. Maybe even more problematically, items 1 and 8 stand apart. Lets try to redo this without these items.


```{r}
reduced_items = item_names[-c(1,8)]
fa.parallel(Q5aar[,reduced_items], cor = "poly")
```

OK, now we only need 3 or maybe only 2 factors. Lets first try with 2.

```{r}
fa_ps = fa(Q5aar[,reduced_items],
nfactors = 2, 
cor = "poly")
fa_ps
```

RMSEA and even more the TLI suggest that 2 factors hould be OK.

This are the new factor loadings: 

```{r}
par(mfrow = c(1,2))
plot(fa_ps)
fa.diagram(fa_ps)
```

This looks like a clear factor structure. Lets look at the items to understand what those factor might measure.

```{r}
factor_loadings = loadings(fa_ps)
simple_loadings = apply(factor_loadings,
1,
function(x)
(abs(x) == max(abs(x))))
item_quest_reduced = item_quest[-c(1,8)]

items_factor1 = which(simple_loadings[1,] == T)
items_factor2 = which(simple_loadings[2,] == T)

cat(paste0("Factor1:\n",
paste(item_quest_reduced[items_factor1],
collapse = "\n"),
"\n"))

cat(paste0("Factor2:\n",
paste(item_quest_reduced[items_factor1],
collapse = "\n"),
"\n"))

```


This looks as if there is one factor around positive reinforcement and communication with the child. Lets call this _positive parenting_. The other factor appears to be about the inability to follow through with threats of punishment, lets call this _inconsistent parenting_ ;-). If we scroll back and look at correlation of the two factors, it looks as if _positive parenting_ has a weak negative correlation with _inconsistent parenting_.


Now, if we want to use the parenting dimensions as predictors, we need to extract the factor scores. Factor scores are likely stored in the `fa_ps` object, which has the results of out exploratory factor analysis. If we want to know what this object contains, we can simply click on in in the "Environment" tab at the right hand side of the RStudio IDE.

```{r}
factor_scores = fa_ps$scores
colnames(factor_scores) = c("positive","inconsistent")
```

Lets quickly look at the factor scores, before we add them to the data.frame with teh 5 year data.:

```{r}
par(mfrow = c(1,2))
hist(factor_scores[,1],
breaks = 25,
main = "positive")
hist(factor_scores[,2],
breaks = 25,
main = "inconsistent")
Q5aar = cbind(Q5aar,factor_scores)

```

If we had a great parenting scale, the histograms should look normally distributed. Especially the histogram for positive parenting deviates from that. It shows that (a) there are a few extreme outliers to the left and (b) there seems to be someting like a ceiling effect. That is, a good number of mothers report to be **VERY** positive in their parenting. This is another instance wehre it would be great to have a social desirability scale in MoBa.

## Some more data wrangling, including calculation of sum scores

Anyhow, lets continue and put together outcome data.
For this, we load the MoBa 8 years data file.

```{r, eval = F}
data_directory = "data/"
Q8aar = read_sav(paste0(data_directory,
"PDB439_Skjema8aar_v10.sav"))
Q8aar = Q8aar[,c("PREG_ID_439","BARN_NR",
paste("NN",c(68:80,        # SMFQ
211:226,      # CCC-2 sort
227:233, 374, # Språåk 20
111:118,      # CD
119:136,      # ADHD
137:144),     # OD
sep = ""))]

Q8aar = data.frame(Q8aar)
```

```{r, echo=F}
load("data/Q8aar.Rdata")
```

Now we are coming to a part, which always takes more time than we would wish: cleaning the data. We wqant to calculate sum-scores, but this is compliated by two facts: MoBa data come with 

* `0` instead of `NA` for data with unambiguous responses
* all scales start at 1
* some of the data are missing.

We will adress these problems in that order, starting with renaming some variables. To make this easier, we'll create simple functions where this is useful.

First we rename all variables:
  
  ```{r}

rename_variables = function(df,old_names,new_names) {
  names(df)[names(df) %in% old_names] = new_names
  return(df)
}

Q8aar = rename_variables(Q8aar,
                         paste0("NN",68:80),
                         paste0("SMFQ.i",1:13,"_8y"))
Q8aar = rename_variables(Q8aar,
                         paste0("NN",c(227:233,374)),
                         paste0("SPRAAK20.i",1:8,"_8y"))
Q8aar = rename_variables(Q8aar,
                         paste0("NN",211:226),
                         paste0("CCCs.i",1:16,"_8y"))
for (item in paste0("CCCs.i",c(10:16),"_8y")) 
  Q8aar[,item] = abs(Q8aar[,item]-5)
Q8aar = rename_variables(Q8aar,
                         paste0("NN",111:118),
                         paste0("CD.i",1:8,"_8y"))
Q8aar = rename_variables(Q8aar,
                         paste0("NN",119:136),
                         paste0("ADHD.i",1:18,"_8y"))
Q8aar = rename_variables(Q8aar,
                         paste0("NN",137:144),
                         paste0("OD.i",1:8,"_8y"))
Q8aar = rename_variables(Q8aar,
                         paste0("NN",145:149),
                         paste0("SCARED.i",1:5,"_8y"))
```

Next we replace all 0 values with NA and subtract 1 from all items. We use the `grep` command, which makes it easy for us to find all variables that match a particular pattern. We also temporarily create a new `data.frame` tmp_items, in which we keep the item data. 

```{r}
all_items = grep("\\.i[0-9]",names(Q8aar), value = T)
tmp_items = Q8aar[,all_items]
tmp_items[tmp_items == 0] = NA
tmp_items = tmp_items-1
```

We can use imputation to replace missing data. I prefer to do this only for cases, for which we have at least 50% of the data. So we remove participants with more than 50% missing data.

```{r}
proportion_missing = rowMeans(is.na(tmp_items))
Q8aar = Q8aar[proportion_missing < .5,]
tmp_items = tmp_items[proportion_missing < .5,]
```

Now we use the function `hotdeck`from the package `VIM` to impute values. (`hotdeck` is fast). Other approaches like nearest neigbourg imputation of multiple chained imputation are more accurate, by they are to slow when analysing bigger data-sets in a turorial.

```{r}
tmp_items_imputed = hotdeck(as.matrix(tmp_items))
tmp_items_imputed = tmp_items_imputed[,1:ncol(tmp_items)]
```

Finally, we can calculate sum scores.

```{r}
scales = c("OD","ADHD","CD","SMFQ","CCC","SPRAAK")
sum_scores = matrix(NA,
                    nrow = nrow(Q8aar),
                    ncol = length(scales))
colnames(sum_scores) = scales
par(mfrow = c(2,3))
for (s in scales) {
  items = grep(s,names(tmp_items_imputed), value = T)
  sum_scores[,s] = rowSums(tmp_items_imputed[,items])
  hist(sum_scores[,s],
       main = s,
       ylab = "")
}
head(sum_scores)
```

This looks more or less as expected. Now we put the Q8aar data together, and merge 5 and 8 years data.

```{r}
Q8aar = cbind(Q8aar, sum_scores)
my_data = merge(Q5aar[,c("PREG_ID_439","BARN_NR",colnames(factor_scores))],
                Q8aar[,c("PREG_ID_439","BARN_NR",scales)],
                by = c("PREG_ID_439","BARN_NR"))
head(my_data)
my_data$positive = as.numeric(scale(my_data$positive))
my_data$inconsistent = as.numeric(scale(my_data$inconsistent))
my_data = my_data[my_data$positive > -4,]
```


## Linear (and other) regressions in R

_For the next section, you can run your own regression by changing the outcome. Use the command_ `names(my_data)` _to see, which variables are available._

Finally, we can do simple linear regression models:
  ```{r}
ADHD_model = lm(ADHD ~ positive + inconsistent, my_data)
summary(ADHD_model)
```

This is exciting! Parenting style at age 5 predicts ADHD symptoms at age 8!
  
  Before getting carried away, it makes sense to look at a few diagnostic plots:
  
  ```{r}
par(ask=F)
plot(ADHD_model)
```

[Here](https://data.library.virginia.edu/diagnostic-plots/) we can go through an explanation of the diagnostic plots.


As expected, this shows that a simple linear regression is probably not appropriate here.

### A beta-binomial regression
A beta-binomial regression would be more appropriate here, because the outcome variable can only be between 0 and an upper bound. There are several possibilities to do such regressions in R. On mature package that can do this is `VGAM`. Alternatives include `brms` and `rstanarm`.

```{r}
library(VGAM)
library(qualityTools)
ADHD_model_bb = vglm(cbind(ADHD, 54-ADHD) ~ poly(positive,2) + poly(inconsistent,2),
                     betabinomial,
                     data = my_data,
                     trace = TRUE)
qqPlot(ADHD_model_bb@fitted.values,"beta", start = list(shape1 = 10, shape2 = 2))
summary(ADHD_model_bb)
```


Even though the QQ plot shows that beta-binomila model is more appropriate, the tutorial continues with the result of the regression model. This is because some of the topics we are going to cover ar easy to implement for a linear regression model, but take a bit more time for alternative models.

### Processing regression results: Tables and plots

With the package `apaTables`we can format results tables from regressions as described in the APA Manual 6.

```{r}
library(apaTables)
apa.reg.table(ADHD_model, filename = "Table2_APA.doc", table.number = 2)
```

It is also relatively easy to plot model results. First, we simply plot the regression coefficients and confidende intervals.

```{r}
library(sjPlot)
library(ggeffects)
library(ggplot2)
plot_model(ADHD_model)

effect_inconsistent = ggpredict(ADHD_model, terms = "inconsistent", typical = "mean")
effect_positive = ggpredict(ADHD_model, terms = "positive", typical = "mean")

ggplot(effect_inconsistent, aes(x = x, y = predicted)) + 
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, x=x, fill = "band"), alpha = 0.3, fill = "blue") + 
  geom_line(col = "blue") + 
  ylab("predicted sum score") + 
  xlab("exposure") + 
  coord_cartesian(ylim = c(0,quantile(my_data$ADHD,c(.95))))
```

We can also plot the effects of teh two predictors together:
  
  ```{r}
my_effects = rbind(cbind(effect_inconsistent,predictor = "inconsistent"),
                   cbind(effect_positive,predictor = "positive"))

ggplot(my_effects, aes(x = x, y = predicted, col = predictor, group = predictor)) + 
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, x=x, fill = predictor, col = NULL), alpha = 0.3) + 
  geom_line() +
  ylab("predicted sum score") + 
  xlab("exposure") + 
  coord_cartesian(ylim = c(0,quantile(my_data$ADHD,c(.95))))
```


### Adding and plotting non-linear effects

This shows a linear effect of parenting style. Is it possible that there are non-linear effects? We can easily test this by adding squared terms to the model.

```{r}
ADHD_model2 = lm(ADHD ~ poly(positive,2,raw = T) + poly(inconsistent,2, raw = T), my_data)
summary(ADHD_model2)
```


We'll use a new library, `cowplot` to put multiple ggplot figures together.

```{r, fig.width= 10}
library(cowplot)
coef_plot = plot_model(ADHD_model2)
effect_inconsistent = ggpredict(ADHD_model2, terms = "inconsistent", typical = "mean")
effect_positive = ggpredict(ADHD_model2, terms = "positive", typical = "mean")

my_effects = rbind(cbind(effect_inconsistent,predictor = "inconsistent"),
cbind(effect_positive,predictor = "positive"))

effect_plot = ggplot(my_effects, aes(x = x, y = predicted, col = predictor, group = predictor)) + 
geom_ribbon(aes(ymin=conf.low, ymax=conf.high, x=x, fill = predictor, col = NULL), alpha = 0.3) + 
geom_line() +
ylab("predicted sum score") + 
xlab("exposure") + 
coord_cartesian(ylim = c(0,quantile(my_data$ADHD,c(.95))))

cowplot::plot_grid(coef_plot,
effect_plot, nrow = 1)
```

### Adding interaction effects

Finally, we can also add an interaction effect:

```{r}
ADHD_model2 = lm(ADHD ~ poly(positive,2,raw = T) + 
poly(inconsistent,2, raw = T) + 
positive:inconsistent,
my_data)
summary(ADHD_model2)

```
This does not seem to make a big difference.

This analysis suggests that mothers who report inconsistent parenting at child age 3 (they do not follow trough with a threatened punishment) later report more ADHD problems with their children. There are several reasons that we should not immediatly endorse a causal interpretation. This could in fact be revers causation, if children with problems at age 8 also had problems at age 5, and it was these problems at age five that caused inconsistent parenting. Some would certainly argue tha this is all genetic: the same genes that are reponsinle for mental health problems also cause inconsistent parenting (e.g. parents with mental health problems might be more prone to be inconsistent.). Genes are an example of confounders: common causes of exposure and outcomes. We do not have genetic information here, but we can easily think of other pontial common causes like parental education and age, parity. In addition, we might want to adjust for other co-variates like the childs gender.


## Investigating potential confounders through plotting

The next lines of R code load these variables and adds them to our data.


```{r, eval = F}
data_directory = "data/"
Q1 = read_sav(paste0(data_directory,
"PDB439_Skjema1_v10.sav"))
Q1 = data.frame(Q1[,c("PREG_ID_439",
paste("AA",c(11,1123:1127,1300:1303),
sep = ""))])
names(Q1)[-1] = c("Q1_YEAR","CivilStatus","mEDUcomp","mEDUcurr",
"fEDUcomp","fEDUcurr","n_people_19p_Q1",
"n_children_12_18_Q1","n_children_6_11_Q1",
"n_children_0_5_Q1"))

Q1$CivilStatus = factor(Q1$CivilStatus,
levels = 1:6,
labels = c("married","separated","cohabitating",
"widow","single","Other"))

EDUlabels = c("basic","1-2 high school","vocational high","general high", "Bachelor","Master" )

Q1$mEDU = Q1$mEDUcomp
idx = is.na(Q1$mEDUcomp) & !is.na(Q1$mEDUcurr) & Q1$mEDUcurr != 1
Q1$mEDU[idx] = Q1$mEDUcurr[idx] - 1
Q1$mEDU = ordered(Q1$mEDU,labels = EDUlabels)
Q1$fEDU = Q1$fEDUcomp
idx = is.na(Q1$fEDUcomp) & !is.na(Q1$fEDUcurr) & Q1$fEDUcurr != 1
Q1$fEDU[idx] = Q1$fEDUcurr[idx] - 1
Q1$fEDU = ordered(Q1$fEDU,labels = EDUlabels)

rm(EDUlabels,idx)


data_directory = "data/"
MFR = read_sav(paste0(data_directory,
"PDB439_MFR_520_v10.sav"))
MFR = data.frame(MFR[,c("PREG_ID_439","BARN_NR","KJONN","MORS_ALDER",
"FAAR","LEVENDEFODTE_5","PARITET_5","FLERFODSEL",
"VEKT","SVLEN_DG","APGAR1","APGAR5","FMND")])
names(MFR)[-c(1,2)] = c("gender","mAge", "birth_year",
"life_births","parity","multiple_births",
"birthweight","pregnancy_duration",
"APGAR1","APGAR5","birthmonth"))
MFR$gender = factor(MFR$gender,levels = 1:2,labels = c("boy","girl"))
MFR[is.na(MFR$life_births), life_births := parity]
MFR[, birthmonth := as.numeric(birthmonth)]


```

```{r, echo = F}
load("data/Q1.Rdata")
load("data/MFR.Rdata")
```

```{r}
my_data = merge(my_data,
MFR,
by = c("PREG_ID_439","BARN_NR"),
all.x = T,
all.y = F)
my_data = merge(my_data,
Q1,
by = c("PREG_ID_439"),
all.x = T,
all.y = F)
```

Now we are readdy to look at a few tables and plots. First, we use the `group_by`function of the `dplyr`package to show the mean exposures grouped by education.

```{r}
library(dplyr)
options(digits = 2)
my_stats = my_data %>% 
group_by(mEDU) %>%
summarise(m = mean(inconsistent),
sd = sd(inconsistent),
N = n()
)
my_stats
```

This does not look very clear. We can also plot the raw data using the nice ggplot package. Lets start with a histogram.

```{r}
library(ggplot2)
ggplot(my_data, aes(x = inconsistent, fill = mEDU)) + 
geom_histogram(aes(y = ..density..),position="dodge2", bins = 10)
```

We can simply plot the `my_stats` table we made above, after adding confidence intervals:


```{r}
my_stats$CIlower = my_stats$m - my_stats$sd/sqrt(my_stats$N)*2
my_stats$CIupper = my_stats$m + my_stats$sd/sqrt(my_stats$N)*2

ggplot(my_stats,aes(x = mEDU, y = m, fill = mEDU)) + geom_bar(stat = "identity") + 
geom_pointrange(aes(x = mEDU, ymin = CIlower, ymax = CIupper)) +
stat_summary(fun.y="median", geom="point", shape=23, size=3, fill="white") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

Now it looks as if the differences between groups are substantial. However, remeber that the expsore variables has an sd of 1, so the differences are maybe not so imporant. To get a better view of the variations within and between groups, we can plot boxplots:

```{r}
ggplot(my_data, aes(x = mEDU, y = inconsistent, fill = mEDU)) + 
geom_boxplot() +
stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

This shows that even though there is are reliable differences between groups, the variation within groups is much larger than the variation between groups.

Now lets look at the association with all potential measured confounders (this looks like a lot of code, but it is really only repetion of the same code!):
```{r, fig.height=10}
breaks = c(seq(15,40, by = 5),50)
my_data$mAgeg = cut(my_data$mAge, breaks = breaks,ordered_result = T)
my_data$parity = as_factor(my_data$parity,ordered = T)
my_data$birthmonth = ordered(my_data$birthmonth)
my_data$birth_year = ordered(my_data$birth_year)

by_edu = ggplot(my_data, aes(x = mEDU, y = inconsistent, fill = mEDU)) + 
geom_boxplot(show.legend = F) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
by_age = ggplot(my_data, aes(x = mAgeg, y = inconsistent, fill = mAgeg)) + 
geom_boxplot(show.legend = F) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
by_parity = ggplot(my_data, aes(x = parity, y = inconsistent, fill = parity)) + 
geom_boxplot(show.legend = F) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
by_gender = ggplot(my_data, aes(x = gender, y = inconsistent, fill = gender)) + 
geom_boxplot(show.legend = F) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
cowplot::plot_grid(by_edu,
by_age,
by_parity,
by_gender)
```

Try to do the same plot but stratify by `birthmonth`.


Here is the same plot for positive parenting:
```{r, fig.height=10}
by_edu = ggplot(my_data, aes(x = mEDU, y = positive, fill = mEDU)) + 
geom_boxplot(show.legend = F) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
by_age = ggplot(my_data, aes(x = mAgeg, y = positive, fill = mAgeg)) + 
geom_boxplot(show.legend = F) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
by_parity = ggplot(my_data, aes(x = parity, y = positive, fill = parity)) + 
geom_boxplot(show.legend = F) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
by_gender = ggplot(my_data, aes(x = gender, y = positive, fill = gender)) + 
geom_boxplot(show.legend = F) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
cowplot::plot_grid(by_edu,
by_age,
by_parity,
by_gender)
```


These associations also do not seem particularily strong.

We conclude this plotting tour with a look at the association between confounders and ADHD symptoms.

```{r}
my_stats = my_data %>% 
group_by(mEDU) %>%
summarise(m = mean(ADHD),
sd = sd(ADHD),
N = n()
)

my_stats$CIlower = my_stats$m - my_stats$sd/sqrt(my_stats$N)*2
my_stats$CIupper = my_stats$m + my_stats$sd/sqrt(my_stats$N)*2

mean_ci_plot = ggplot(my_stats,aes(x = mEDU, y = m, fill = mEDU)) + 
geom_bar(stat = "identity", show.legend = F) + 
geom_pointrange(aes(x = mEDU, ymin = CIlower, ymax = CIupper), show.legend = F)  + 
ylab("ADHD") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

box_plot = ggplot(my_data, aes(x = mEDU, y = ADHD, fill = mEDU)) + 
geom_boxplot(show.legend = F) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

cowplot::plot_grid(mean_ci_plot,
box_plot)
```

OK, now that we have seen that exposure and outcome are (weakly) associated with a likely confounder, it makes sense to include them into the model. 

## A final regression model: Ordinal predictors

```{r}
ADHD_model4 = lm(ADHD ~
poly(positive,2, raw = T) + 
poly(inconsistent,2, raw = T) + 
mEDU + 
poly(mAge,2, raw = T) + 
parity +
gender,
data = my_data)
summary(ADHD_model4)
```

This looks a bit wild. The reason is that if we include ordinal variables like education and parity into the model, R will by default implement polynomials to a degree equal to `number_levels - 1`. However, we can correct this by setting the contrast to for example a polynomial of order 2.

```{r}
n_levels_Edu = length(levels(my_data$mEDU))
n_levels_parity = length(levels(my_data$parity))
n_levels_bm = length(levels(my_data$birthmonth))
contrasts(my_data$mEDU, 2) <- contr.poly(n_levels_Edu)
contrasts(my_data$parity, 2) <- contr.poly(n_levels_parity)
contrasts(my_data$birthmonth, 2) <- contr.poly(n_levels_bm)

ADHD_model4 = lm(ADHD ~
poly(positive,2, raw = T) + 
poly(inconsistent,2, raw = T) + 
mEDU + 
poly(mAge,2, raw = T) + 
parity +
gender + 
birthmonth,
data = my_data)
summary(ADHD_model4)
plot_model(ADHD_model4)
```

This table suggests that parenting style still has an effect, after we controlled for a number of potential confounders. It is, however, difficult to understand the importance (effect sizes) for the different predictors, because they are measured on different scales (have different means and standard deviations). To get a better idea of effect sizesm we can again plot marginal effects.

We use a small loop, so that we don't have to write the same lines of code repeatedly for the different predictors.

```{r}
ylim = c(0,quantile(my_data$ADHD,c(.95)))

predictors = c("positive","inconsistent","mAge","mEDU","parity","gender", "birthmonth")
effect_plots = vector(mode = "list", length = length(predictors))
names(effect_plots) = predictors

for ( p in 1:length(predictors)) {
  efct = ggpredict(ADHD_model4, terms = predictors[p])
  effect_plots[[p]] = ggplot(efct, aes(x = x, y = predicted)) + 
    geom_ribbon(aes(ymin=conf.low,
                    ymax=conf.high,
                    x=x,
                    col = NULL),
                alpha = 0.3) + 
    geom_line() +
    ylab("predicted sum score") + 
    xlab(predictors[p]) + 
    coord_cartesian(ylim = ylim)
}
```

First we look at the effect of confounders:
  
  ```{r}
cowplot::plot_grid(plotlist = effect_plots[3:length(predictors)])
```

And finally on the effect of our exposures:
  
  ```{r}
cowplot::plot_grid(plotlist = effect_plots[1:2])

```

Interestingly, this suggests that the association between inconsisten parenting style at child age 5 and ADHD sum score at child age 8 increases after we control for some potential confounders.

library(ggExtra)
p = ggplot(my_data, aes(x=inconsistent, y=ADHD)) +
  geom_hex(stat = "binhex") + 
  scale_fill_gradientn(colours=c("white","black"),name = "Frequency") + 
  geom_smooth(method = "lm") +
  coord_cartesian(ylim = c(0,quantile(my_data$ADHD,c(.975))))


p = ggMarginal(p, type="histogram")
p
