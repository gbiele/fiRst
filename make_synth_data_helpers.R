library(GenOrd)
source("ordsample.R")
source("ordcont.R")
source("corrcheck.R")
library(stringr)
library(VIM)
library(haven)
library(foreign)

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
  for (v in grep("PREG|BARN",names(df), invert = T, value = T)) {
    x = gsub(paste0('\"',v,'\"'),paste('\"',attr(df[,v],"label"),'\"'),x)
  }
  
  cat(x,file = fn,sep = "\n")
}



write2spss = function(df,name, fix_variable_labels = T) {
  fn_csv = paste0(name,".csv")
  fn_sps = paste0(name,".sps")
  
  foreign:::writeForeignSPSS(df,
                             datafile = fn_csv,
                             codefile = fn_sps)
  
  if (fix_variable_labels) {
    tmp = read.csv(fn_csv)
    tmp = tmp[,1:ncol(df)]
    names(tmp) = names(df)[1:ncol(tmp)]
    for (v in grep("PREG|BARN",names(df), invert = T, value = T)) {
      if (any(class(df[,v]) == "factor")) {
        if (length(grep("kryss",levels(df[,v])[1])) > 0)
          tmp[,v] = tmp[,v]-1
      }
              
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
}

clean_vars = function(df) {
  df = data.frame(df)
  for (v in names(df)) {
    if (class(df[,v]) %in% c("labelled","haven_labelled")) {
      lbls = attr(df[,v],"labels")
      is_multicross_val = lbls["Mer enn ett kryss"]
      if ( length(is_multicross_val) > 0 ) {
        df[which(df[,v] == is_multicross_val),v] = NA
      } 
    }
  }
  return(df)
}

make_synth_data = function(df, N = NULL) {
  
  if (is.null(N)) N = nrow(df)
  
  margins = lapply(df, function(x) {
    p = prop.table(table(x))
    if (min(p) < 1e-4) {
      delta = 1e-4-p[p<1e-4]
      p[p<1e-4] = p[p<1e-4]+delta
      p[p>1e-3] = p[p>1e-3]-sum(delta)/length(p[p>1e-3])
    }
    m = cumsum(p)
    return(m[-length(m)])
  })
  
  Sigma = cor(df, use = "pairwise.complete.obs",method = "spearman")
  diag(Sigma) = 1
  
  # synthesize
  sdf = ordsample(nrow(df),
                  marginal = margins,
                  Sigma = Sigma,
                  Spearman = T)
  
  sdf = data.frame(sdf)
  
  for (v in names(df))
    sdf[is.na(df[,v]),v] = NA
  
  sdf = sdf[sample(nrow(df),N),]
  
  smargins = lapply(sdf, function(x) {
    m = cumsum(prop.table(table(x)))
    return(m[-length(m)])
  })
  
  ## correct for slightly off margins
  off_margins = names(which(sapply(margins,length) != sapply(smargins, length)))
  for (v in off_margins) {
    m = table(df[,v])
    sm = table(sdf[,v])
    fs = merge(data.frame(m),
               data.frame(sm),
               by = "Var1",
               all = T)
    fs$Var1 = as.numeric(fs$Var1)
    missing_val = fs[is.na(fs$Freq.y),"Var1"]
    nearest_val = fs$Var1[abs((fs$Var1-missing_val)) == 1]
    switch_idx = sample(which(sdf[,v] == nearest_val),1)
    sdf[switch_idx,v] = missing_val
  }
  
  # statistics for symthesized data
  sSigma = cor(sdf, use = "pairwise.complete.obs",method = "spearman")
  diag(sSigma) = 1
  smargins = lapply(sdf, function(x) {
    m = cumsum(prop.table(table(x)))
    return(m[-length(m)])
  })
  
  par(mfrow = c(2,2),mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
  plot(Sigma,sSigma)
  hist(Sigma-sSigma, main = "delta Sigma")
  plot(do.call(c,margins),do.call(c,smargins))
  hist(do.call(c,margins)-do.call(c,smargins), main = "delta margins")
  return(sdf)
}

add_label_s = function(df,labelled_df) {
  for (v in names(labelled_df)) {
    if(class(labelled_df[,v])  %in% c("labelled","haven_labelled")) {
      df[,v] = labelled(df[,v],labels = attr(labelled_df[,v],"labels"))
      df[,v] = as_factor(df[,v], ordered = T)
      attr(df[,v],"label") = attr(labelled_df[,v],"label")
    }
  }
  return(df)
}