source("make_synth_data_helpers.R")

data_directory = "F:/Forskningsprosjekter/PDB 439 - Language and learnin_/Forskningsfiler/1-OrginalDATA-kunKOPIEREfra/1-Pdb439_MoBa_v10/Data/"


Q1 = read_sav(paste0(data_directory,
                     "PDB439_Skjema1_v10.sav"))
Q1 = data.frame(Q1[,c("PREG_ID_439",
                      paste0("AA",c(1124:1125)))])
Q1 = clean_vars(Q1)

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
Q5aar = Q5aar[rowSums(is.na(Q5aar[,grep("LL",names(Q5aar))])) < 5,]
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
ALL = merge(ALL,Q1, by = c("PREG_ID_439"))

sMFRQ1 = merge(MFR,Q1, by = c("PREG_ID_439"))

sALL = make_synth_data(ALL[,-c(1,2)])
sQ5aar = make_synth_data(Q5aar[,-c(1,2)],N = nrow(Q5aar)-nrow(ALL))
#sQ8aar = make_synth_data(Q8aar[,-c(1,2)],N = nrow(Q8aar)-nrow(ALL))
sMFR = make_synth_data(MFR[,-c(1,2)],N = nrow(MFR)-nrow(ALL)-nrow(Q5aar)-nrow(Q8aar))
sMFRQ1 = make_synth_data(sMFRQ1[,-c(1,2)],N = nrow(sMFRQ1)-nrow(ALL)-nrow(Q5aar)-nrow(Q8aar))


for (v in setdiff(names(sALL),names(sQ5aar)))
  sQ5aar[,v] = NA
for (v in setdiff(names(sALL),names(sMFR)))
  sMFR[,v] = NA

sALL = rbind(sALL,sMFR)
sALL = rbind(sALL,sQ5aar)
#sALL = rbind(sALL,sQ8aar)

sALL = add_label_s(sALL,Q8aar)
sALL = add_label_s(sALL,Q5aar[,intersect(names(Q5aar),names(sALL))])
sALL = add_label_s(sALL,MFR[,intersect(names(MFR),names(sALL))])

sALL$PREG_ID_439 = 1:nrow(sALL)
sALL$BARN_NR = sample(sort(unique(ALL$BARN_NR)),
                      nrow(sALL),
                      replace = T,
                      prob = prop.table(table(ALL$BARN_NR)))

q5 = sALL[,grep("PREGID|BARN_NR|^LL",names(sALL))]
q8 = sALL[,grep("PREGID|BARN_NR|^NN",names(sALL))]
mfr = sALL[,grep("^LL|^NN",names(sALL), invert = T)]

write2spss(q5,"H:/misc/fiRst/data/q5")
write2spss(q8,"H:/misc/fiRst/data/q8")
write2spss(mfr,"H:/misc/fiRst/data/mfr")

