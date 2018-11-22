library(GenOrd)
library(VIM)

load("data/Q5aar.Rdata")
Q5aar_orig = Q5aar

item_names = names(Q5aar)[-c(1,2)]

for (i in item_names)
  Q5aar[which(Q5aar[,i] == 0),i] = NA

is_missing = is.na(Q5aar[,item_names])
Q5aar = Q5aar[rowMeans(is_missing[,item_names]) < .5,]
is_missing = is.na(Q5aar[,item_names])
ps_imputations = hotdeck(as.matrix(Q5aar[,item_names]))
Q5aar[,item_names] = ps_imputations[,item_names]

is_missing = is_missing[rowSums(Q5aar[,item_names]) > 9,]
Q5aar = Q5aar[rowSums(Q5aar[,item_names]) > 9,]

n = nrow(Q5aar)
margins = lapply(Q5aar[,3:11], function(x) {
  return(head(cumsum(prop.table(table(x))),4))
})
Sigma = cor(Q5aar[,item_names], use = "pairwise.complete.obs")

x = ordsample(round(n*1.1),
          marginal = margins,
          Sigma = Sigma)

x = cbind(x,
          F1 = rowMeans(x[,c(2,4,9)]),
          F2 = rowMeans(x[,c(3,5,6,7)]))

Q5aar = cbind(Q5aar,
              F1 = rowMeans(Q5aar[,item_names[c(2,4,9)]]),
              F2 = rowMeans(Q5aar[,item_names[c(3,5,6,7)]]))

my_cov = cov(Q5aar[,c("F1","F2")])
xF = as.matrix(x[,c("F1","F2")])
synth_ps = matrix(NA,
                  nrow = nrow(Q5aar),
                  ncol = 9)
tmp_x = x
for (k in 1:nrow(Q5aar)) {
  idx = which.min(mahala(xF,
                            as.matrix(Q5aar[k,c("F1","F2")]),
                            sigma = my_cov))
  synth_ps[k,] = as.matrix(tmp_x[idx,1:9])
  tmp_x = tmp_x[-idx,]
  xF = xF[-idx,]
}

synth_ps = cbind(synth_ps,
          F1 = rowMeans(synth_ps[,c(2,4,9)]),
          F2 = rowMeans(synth_ps[,c(3,5,6,7)]))

smoothScatter(synth_ps[,"F1"],Q5aar[,"F1"])
smoothScatter(synth_ps[,"F2"],Q5aar[,"F2"])

colnames(synth_ps)[1:9] = item_names

for (v in item_names) {
  Q5aar[,v] = synth_ps[,v]
  Q5aar[is_missing[,v],v] = NA
  Q5aar[,v] = labelled(Q5aar[,v],attr(Q5aar_orig[,v], "labels"))
  attr(Q5aar[,v], "label") = attr(Q5aar_orig[,v], "label")
}
Q5aar = Q5aar[,1:11]
save(Q5aar,file = "data/sQ5aar.Rdata")


load("data/Q8aar.Rdata")

scale_items = list(SMFQ = 68:80,
                   CCC = 211:226,
                   SPRAAK20 = c(227:233, 374),
                   CD = 111:118,
                   ADHD = 119:136,
                   OD = 137:144)

for (s in names(scale_items)) {
  items = paste0("NN",scale_items[[s]])
  for (k in 1:nrow(Q8aar)) {
    Q8aar[k,items] = sample(Q8aar[k,items])
  }
}

save(Q8aar, file = "data/sQ8aar.Rdata")

load("data/MFR.Rdata")
MFR$mAge = MFR$mAge + sample(-2:2,nrow(MFR),replace = T)
save(MFR,file = "data/sMFR.Rdata")
