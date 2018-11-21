library(GenOrd)

load("data/Q5aar.Rdata")
items = grep("LL",names(Q5aar))
for (i in items)
  Q5aar[which(Q5aar[,i] == 0),i] = NA
Q5aar = Q5aar[complete.cases(Q5aar),]
Q5aar = Q5aar[rowSums(Q5aar[,items]) > 9,]

n = nrow(Q5aar)
margins = lapply(Q5aar[,3:11], function(x) {
  return(head(cumsum(prop.table(table(x))),4))
})
Sigma = cor(Q5aar[,items], use = "pairwise.complete.obs")
margins[5]$LL530 = margins[5]$LL530

x = ordsample(n,
          marginal = margins,
          Sigma = Sigma)

x = cbind(x,
          F1 = rowMeans(x[,c(2,4,9)]),
          F2 = rowMeans(x[,c(3,5,6,7)]))

Q5aar = cbind(Q5aar,
              F1 = rowMeans(Q5aar[,items[c(2,4,9)]]),
              F2 = rowMeans(Q5aar[,items[c(3,5,6,7)]]))

x = data.table(x)
Q5aar = data.table(Q5aar)

setkeyv(x,c("F2","F1"))
setkeyv(Q5aar,c("F2","F1"))

smoothScatter(x$F1,Q5aar$F1)
smoothScatter(x$F2,Q5aar$F2)

cor(x$LL527,Q5aar$LL527)
smoothScatter(x$F2,Q5aar$F2)