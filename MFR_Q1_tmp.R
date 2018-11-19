library(haven)
library(data.table)
library(mvnfast)
MoBaroot = "F:/Forskningsprosjekter/PDB 439 - Language and learnin_/Forskningsfiler/1-OrginalDATA-kunKOPIEREfra/1-Pdb439_MoBa_v10/Data/"


MFR = read_sav(paste(MoBaroot,"PDB439_MFR_520_v10.sav",sep = ""))
MFR = data.table(MFR[,c("PREG_ID_439","BARN_NR","KJONN","MORS_ALDER",
                        "FAAR","LEVENDEFODTE_5","PARITET_5","FLERFODSEL",
                        "VEKT","SVLEN_DG","APGAR1","APGAR5","FMND")])
setnames(MFR,names(MFR)[-c(1,2)],c("gender","mAge",
                                   "birth_year","life_births","parity","multiple_births",
                                   "birthweight","pregnancy_duration","APGAR1","APGAR5","birthmonth"))
MFR$gender = factor(MFR$gender,levels = 1:2,labels = c("boy","girl"))
MFR[is.na(life_births), life_births := parity]
MFR[, birthmonth := as.numeric(birthmonth)]

save(MFR,file = "data/MFR.Rdata")


Q1 = read_sav(paste(MoBaroot,"PDB439_Skjema1_v10.sav",sep = ""))
Q1 = data.table(Q1[,c("PREG_ID_439",paste("AA",c(11,1123:1127,1300:1303),sep = ""))])
Q1 = clean_vars(Q1)
setnames(Q1,names(Q1)[-1],c("Q1_YEAR","CivilStatus","mEDUcomp","mEDUcurr","fEDUcomp","fEDUcurr","n_people_19p_Q1","n_children_12_18_Q1","n_children_6_11_Q1","n_children_0_5_Q1"))
Q1$CivilStatus = factor(Q1$CivilStatus,levels = 1:6,labels = c("married","separated","cohabitating","widow","single","Other"))

EDUlabels = c("basic","1-2 high school","vocational high","general high", "Bachelor","Master" )

Q1$mEDU = Q1$mEDUcomp
idx = is.na(Q1$mEDUcomp) & !is.na(Q1$mEDUcurr) & Q1$mEDUcurr != 1
Q1$mEDU[idx] = Q1$mEDUcurr[idx] - 1
Q1$mEDU = ordered(Q1$mEDU,labels = EDUlabels)
Q1$fEDU = Q1$fEDUcomp
idx = is.na(Q1$fEDUcomp) & !is.na(Q1$fEDUcurr) & Q1$fEDUcurr != 1
Q1$fEDU[idx] = Q1$fEDUcurr[idx] - 1
Q1$fEDU = ordered(Q1$fEDU,labels = EDUlabels)

Q1[is.na(n_children_0_5_Q1) & (!is.na(n_children_6_11_Q1) | !is.na(n_children_12_18_Q1 | !is.na(n_people_19p_Q1))),n_children_0_5_Q1 := 0]
Q1[is.na(n_children_6_11_Q1) & (!is.na(n_children_0_5_Q1) | !is.na(n_children_12_18_Q1 | !is.na(n_people_19p_Q1))),n_children_6_11_Q1 := 0]
Q1[is.na(n_children_12_18_Q1) & (!is.na(n_children_6_11_Q1) | !is.na(n_children_0_5_Q1 | !is.na(n_people_19p_Q1))),n_children_12_18_Q1 := 0]
Q1[,n_children_Q1 := n_children_0_5_Q1 + n_children_6_11_Q1 + n_children_12_18_Q1]
Q1[n_children_Q1 > 10,n_children_Q1 := 10]
Q1[,Q1 := T]

rm(EDUlabels,idx)

save(Q1,file = "data/Q1.Rdata")