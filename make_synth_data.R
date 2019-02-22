library(GenOrd)
library(VIM)
library(haven)
library(foreign)

#################################################
################ helper functions ###############
#################################################

set_value_labels = function(fn) {
  y = readLines(fn)
  valabidx = grep("VALUE LABELS",y):(grep("VARIABLE LEVEL",y)[1]-1)
  y[valabidx] = gsub('1 \"','0 \"',y[valabidx])
  y[valabidx] = gsub('2 \"','1 \"',y[valabidx])
  y[valabidx] = gsub('3 \"','2 \"',y[valabidx])
  y[valabidx] = gsub('4 \"','3 \"',y[valabidx])
  y[valabidx] = gsub('5 \"','4 \"',y[valabidx])
  y[valabidx] = gsub('6 \"','5 \"',y[valabidx])
  cat(y,file = fn,sep = "\n")
}

set_variable_labels = function(df,fn) {
  x = readLines(fn)
  
  for (v in names(df)[-c(1,2)])
    x = gsub(paste0('"',v,'"'),
             paste0('"',attr(df[,v],"label"),'"'),
             x)
  
  cat(x,file = fn,sep = "\n")
}


set_save_file = function(df,fn) {
  x = readLines(fn)
  end_idx = grep("EXECUTE",x)
  save_string = paste0('SAVE OUTFILE = \"',gsub('sps','sav',fn),'\"')
  x = c(x[end_idx-1],
        save_string,
        x[end_idx:length(x)])
  cat(x,file = fn,sep = "\n")
}

set_variable_labels = function(df,fn) {
  x = readLines(fn)
  
  for (v in names(df)[-c(1,2)])
    x = gsub(paste0('"',v,'"'),
             paste0('"',attr(df[,v],"label"),'"'),
             x)
  
  cat(x,file = fn,sep = "\n")
}


write2spss = function(df,name, fix_variable_labels = T) {
  fn_csv = paste0(name,".csv")
  fn_sps = paste0(name,".sps")
  
  for (v in names(df)) {
    lbl = attr(df[,v],"label")
    if (!is.null(lbl)) {
      lbls = names(attr(df[,v],"labels"))
      if (!is.null(lbls)) {
        if (length(grep("kryss",lbls[1])) > 0) {
          df[,v] = ordered(df[,v],
                           levels = 0:(length(lbls)-1),
                           labels = lbls)
          attr(df[,v],"label") = lbl
        }
      }
    }
  }
    
  
  foreign:::writeForeignSPSS(df,
                             datafile = fn_csv,
                             codefile = fn_sps)
  
  if (fix_variable_labels) {
    tmp = read.csv(fn_csv)
    names(tmp)[1:length(names(df))] = names(df)
    for (v in names(which(sapply(df, is.factor)))) {
      if (length(grep("kryss",levels(df[,v])[1])) > 0)
         tmp[,v] = tmp[,v]-1  
    }
    
    
    write.table(tmp,
                file = fn_csv,
                na = "",
                col.names = F,
                row.names = F,
                sep = ",")
    set_variable_labels(df,fn_sps)
  }
  
  
  set_value_labels(fn_sps)
  
  set_save_file(fn_sps)
}

#################################################
##################### Q5aar #####################
#################################################


load("data/real_data/Q5aar.Rdata")
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

tmp = which(rowSums(Q5aar_orig == 0, na.rm = T) > 0)
for (k in tmp) {
  o = Q5aar_orig[k,]
  if (sum(Q5aar$PREG_ID_439 == o$PREG_ID_439 & Q5aar$BARN_NR == o$BARN_NR) > 0) {
    idx = which(Q5aar$PREG_ID_439 == o$PREG_ID_439 & Q5aar$BARN_NR == o$BARN_NR)
    Q5aar[idx,which(o == 0)] = 0
  }
}

save(Q5aar,file = "data/sQ5aar.Rdata")

write2spss(Q5aar,"H:/misc/fiRst/data/sQ5aar")

#################################################
##################### Q8aar #####################
#################################################

load("data/real_data/Q8aar.Rdata")

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
write2spss(Q8aar,"data/sQ8aar")

#################################################
###################### MFR ######################
#################################################

clean_vars = function(df) {
  df = data.frame(df)
  for (v in names(df)) {
    if (class(df[,v]) == "labelled") {
      lbls = attr(df[,v],"labels")
      is_multicross_val = lbls["Mer enn ett kryss"]
      if ( length(is_multicross_val) > 0 ) {
        df[which(df[,v] == is_multicross_val),v] = NA
      } 
    }
  }
  return(df)
}

data_directory = "F:/Forskningsprosjekter/PDB 439 - Language and learnin_/Forskningsfiler/1-OrginalDATA-kunKOPIEREfra/1-Pdb439_MoBa_v10/Data/"

MFR = read_sav(paste0(data_directory,
                        "PDB439_MFR_520_v10.sav"))
MFR = data.frame(MFR[,c("PREG_ID_439","BARN_NR",
                        "KJONN","MORS_ALDER",
                        "PARITET_5","FMND")])
ma = cut(MFR$MORS_ALDER, breaks = seq(10,50,by = 5))
lbs = 1:max(as.numeric(ma))
names(lbs) = levels(ma)
MFR$MORS_ALDER = labelled(x = as.numeric(ma), labels = lbs)
MFR$FMND = as.numeric(MFR$FMND)

Q5aar = read_sav(paste0(data_directory,
                        "PDB439_Skjema5aar_v10.sav"))
Q5aar = Q5aar[,c("PREG_ID_439","BARN_NR",
                 paste0("LL",526:534))]
Q5aar = clean_vars(Q5aar)

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
Q8aar = clean_vars(Q8aar)

ALL = merge(Q5aar,MFR, by = c("PREG_ID_439","BARN_NR"))
ALL = merge(ALL,Q8aar, by = c("PREG_ID_439","BARN_NR"))
ALL = ALL[rowSums(is.na(ALL[,grep("LL",names(ALL))])) < 5,]

ALLd = ALL[,-c(1,2)]


margins = lapply(ALLd, function(x) {
  m = cumsum(prop.table(table(x)))
  return(m[m != 1])
})

Sigma = cor(ALLd, use = "pairwise.complete.obs",method = "spearman")
diag(Sigma) = 1

k = 37
x = ordsample(1000,
              marginal = margins[1:k],
              Sigma = Sigma[1:k,1:k],
              Spearman = T)
