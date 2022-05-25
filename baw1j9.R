#install have package if not installed 
#package(haven)
library(haven)

dir()
 #check working directory

adva<-read_xpt("c4591001-A-D-adva.zip")
adva<-as.data.frame(adva)
 #check size
dim(adva)
 #check number of subjects
length(unique(adva$SUBJID))


table<-matrix(c(0,0,0,0),nrow=2,ncol=2,dimnames = list(c("POS","NEG"),c("Placebo","BNT162b2")))

nabdata_2020<-nabdata[nabdata$ISDTC<"2020-11-15",]
length(unique(nabdata_2020$SUBJID))
nabdata_2020_visit3<-nabdata_2020[nabdata_2020$VISIT=="V3_MONTH1_POSTVAX2_L",]
#cbind(nabdata_2020_visit3$SUBJID[1:100], nabdata_2020_visit3$ISDTC[1:100], nabdata_2020_visit3$VISIT[1:100], nabdata_2020_visit3$TRTA[1:100])
length(unique(nabdata_2020_visit3$SUBJID))
table[1,]<-c(sum(nabdata_2020_visit3$TRTA=="Placebo"),sum(nabdata_2020_visit3$TRTA=="BNT162b2 Phase 2/3 (30 mcg)"))

#write the table to a csv file ---------------------------------
#write.csv(nabdata_2020_visit3, "visit3_NAB_positive.csv")


nabdata<-adva[(adva$PARAM=="N-binding antibody - N-binding Antibody Assay" & adva$AVALC=="NEG"),]
nabdata_2020<-nabdata[nabdata$ISDTC<"2020-11-15",]
length(unique(nabdata_2020$SUBJID))
nabdata_2020_visit3<-nabdata_2020[nabdata_2020$VISIT=="V3_MONTH1_POSTVAX2_L",]
length(unique(nabdata_2020_visit3$SUBJID))
table[2,]<-c(sum(nabdata_2020_visit3$TRTA=="Placebo"),sum(nabdata_2020_visit3$TRTA=="BNT162b2 Phase 2/3 (30 mcg)"))

table

#============================================

#repeat for initial visit


table<-matrix(c(0,0,0,0),nrow=2,ncol=2,dimnames = list(c("POS","NEG"),c("Placebo","BNT162b2")))

nabdata<-adva[(adva$PARAM=="N-binding antibody - N-binding Antibody Assay" & adva$AVALC=="POS"),]
nabdata_2020<-nabdata[nabdata$ISDTC<"2020-11-15",]
length(unique(nabdata_2020$SUBJID))
nabdata_2020_visit1<-nabdata_2020[nabdata_2020$VISIT=="V1_DAY1_VAX1_L",]
#cbind(nabdata_2020_visit1$SUBJID[1:100], nabdata_2020_visit1$ISDTC[1:100], nabdata_2020_visit1$VISIT[1:100], #nabdata_2020_visit1$TRTA[1:100])
length(unique(nabdata_2020_visit1$SUBJID))
table[1,]<-c(sum(nabdata_2020_visit1$TRTA=="Placebo"),sum(nabdata_2020_visit1$TRTA=="BNT162b2 Phase 2/3 (30 mcg)"))

nabdata<-adva[(adva$PARAM=="N-binding antibody - N-binding Antibody Assay" & adva$AVALC=="NEG"),]
nabdata_2020<-nabdata[nabdata$ISDTC<"2020-11-15",]
length(unique(nabdata_2020$SUBJID))
nabdata_2020_visit1<-nabdata_2020[nabdata_2020$VISIT=="V1_DAY1_VAX1_L",]
length(unique(nabdata_2020_visit1$SUBJID))
table[2,]<-c(sum(nabdata_2020_visit1$TRTA=="Placebo"),sum(nabdata_2020_visit1$TRTA=="BNT162b2 Phase 2/3 (30 mcg)"))

table


#==============================
#create table with visits cross tabulated

nabdata<-adva[(adva$PARAM=="N-binding antibody - N-binding Antibody Assay"),]
nabdata_2020<-nabdata[nabdata$ISDTC<"2020-11-15",]
length(unique(nabdata_2020$SUBJID))
nabdata_2020_visit1<-nabdata_2020[nabdata_2020$VISIT=="V1_DAY1_VAX1_L",]
nabdata_2020_visit3<-nabdata_2020[nabdata_2020$VISIT=="V3_MONTH1_POSTVAX2_L",]

  subjects<-unique(nabdata_2020_visit$SUBJID)

visit1_match<-match(subjects,nabdata_2020_visit1$SUBJID)
visit3_match<-match(subjects,nabdata_2020_visit3$SUBJID)


nabtable<- as.data.frame(cbind(subjects, nabdata_2020_visit1$SUBJID[visit1_match], nabdata_2020_visit1$TRTA[visit1_match], nabdata_2020_visit1$ISDTC[visit1_match], nabdata_2020_visit1$VISIT[visit1_match], nabdata_2020_visit1$PARAM[visit1_match], nabdata_2020_visit1$AVALC[visit1_match], nabdata_2020_visit3$SUBJID[visit3_match], nabdata_2020_visit3$TRTA[visit3_match], nabdata_2020_visit3$ISDTC[visit3_match], nabdata_2020_visit3$VISIT[visit3_match], nabdata_2020_visit3$PARAM[visit3_match], nabdata_2020_visit3$AVALC[visit3_match])
colnames(nabtable)<-c("SUBJID","SUBJID_V1","TRTA_V1","ISDTC_V1","VISIT_V1","PARAM_V1","AVALC_V1", "SUBJID_V3","TRTA_V3","ISDTC_V3","VISIT_V3","PARAM_V3","AVALC_V3")


#display results

sum((nabtable$AVALC_V1=="NEG" & nabtable$TRTA_V1=="Placebo" & nabtable$AVALC_V3=="POS"), na.rm = TRUE)
sum((nabtable$AVALC_V1=="NEG" & nabtable$TRTA_V1=="Placebo" & nabtable$AVALC_V3=="NEG"), na.rm = TRUE)
sum((nabtable$AVALC_V1=="NEG" & nabtable$TRTA_V1=="BNT162b2 Phase 2/3 (30 mcg)" & nabtable$AVALC_V3=="POS"), na.rm = TRUE)
sum((nabtable$AVALC_V1=="POS" & nabtable$TRTA_V1=="BNT162b2 Phase 2/3 (30 mcg)" & nabtable$AVALC_V3=="NEG"), na.rm = TRUE)
sum((nabtable$AVALC_V1=="POS" & nabtable$TRTA_V1=="Placebo" & nabtable$AVALC_V3=="NEG"), na.rm = TRUE)
sum((nabtable$AVALC_V1=="POS" & nabtable$TRTA_V1=="BNT162b2 Phase 2/3 (30 mcg)" & nabtable$AVALC_V3=="POS"), na.rm = TRUE)
sum((nabtable$AVALC_V1=="POS" & nabtable$TRTA_V1=="Placebo" & nabtable$AVALC_V3=="POS"), na.rm = TRUE)
sum((nabtable$AVALC_V1=="NEG" & nabtable$TRTA_V1=="BNT162b2 Phase 2/3 (30 mcg)" & nabtable$AVALC_V3=="NEG"), na.rm = TRUE)
sum((nabtable$AVALC_V1=="NEG" & nabtable$TRTA_V1=="Placebo" & nabtable$AVALC_V3=="NEG"), na.rm = TRUE)



