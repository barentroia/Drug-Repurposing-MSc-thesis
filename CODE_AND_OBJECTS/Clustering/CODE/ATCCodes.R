##Making of the ATC codes and DrugTargets dataframe - these dataframes match CMap drugs to ATC codes/MoA
#Drug data was retrieved from ChEMBL and from the atc package
library(reshape2)
library(VennDiagram)

ChemblData=read.csv("../DATA/atccodes.csv",header = TRUE)
load("../DATA/AllData-WHOCC-dump-2016-02-12.RData")
load("../../Replication_DN_Iorioetal2010/DATA/DRUG_DISTANCES.ro")
DrugTargets=read.csv("../DATA/drugtargestschembl.csv",header=TRUE)

#Step1 - Only retrieve the information in ChEMBL data relative to the drugs in Cmap data
#Match the drug identifiers in DRUG_DISTANCES to one of the synonyms in Chembl data. 
match=as.vector(sapply(colnames(DRUG_DISTANCES),function(x){
  string= paste("^(.*; )?\\b",x,"\\b( \\(.*)?$",sep="")
  matches=grep(string,ChemblData$SYNONYMS,ignore.case=TRUE) 
}))
ChData_cmap=melt(match)
colnames(ChData_cmap)=c("row_orig_ChData","Cmap_drugs")

duplic=unique(ChData_cmap$Cmap_drugs[duplicated(ChData_cmap$Cmap_drugs)]) #only 2 drugs, had duplicated entries in the original data
ChData_cmap=ChData_cmap[-c(333,497),]

ATC_cmap=cbind(ChData_cmap,ChemblData[ChData_cmap$row_orig_ChData,])
ATC_cmap=ATC_cmap[which(ATC_cmap$ATC_CODE!=""),]  #This dataframe is the reference data frame.

#From this table, we can retrieve only the correspondence atc code to drug (where one drug can have multiple atc codes)
atcs <- strsplit(as.character(ATC_cmap$ATC_CODE), '; ')
ATCcodes=data.frame(Drugs=rep(ATC_cmap$Cmap_drugs,sapply(atcs, FUN=length)),ATC_codes=unlist(atcs))

#Step 2 - Retrieve data from WHOCC for cmap drugs
atc <- as.data.frame(AllData[["atc"]], stringsAsFactors=F)
atc <- atc[nchar(atc$key) == 7, ]
atc <- atc[atc$name %in% colnames(DRUG_DISTANCES), ] #731 drugs

#Step 3 - See if the common codes correspond to the same drugs in the two datasets
Drugs_ch=as.character(ATCcodes$Drugs[ATCcodes$ATC_codes%in%atc$key])
Drugs_atc=atc$name[match(ATCcodes$ATC_codes,atc$key)]
Drugs_atc=Drugs_atc[!is.na(Drugs_atc)]

drugs=list(ChEMBL=ATCcodes$Drugs,WHOCC=atc$name)
venn.diagram(drugs,filename="Venn_drugs.tiff")

#I am now adding the info related to the drugs than only appear in ChEMBL data:
chDrugs=ATCcodes$Drugs[!ATCcodes$Drugs%in%atc$name]
ATCcodes[ATCcodes$Drugs%in%as.character(chDrugs),]

colnames(atc)=c("Codes","Drugs")
ATCcodes=ATCcodes[,c(2,1)]
colnames(ATCcodes)=c("Codes","Drugs")
MergedATCcodes=rbind(atc,ATCcodes[ATCcodes$Drugs%in%as.character(chDrugs),])

write.csv(MergedATCcodes,"../DATA/ProcessedATCcodes.csv",row.names = FALSE)

#Finding the MoA of Cmap drugs (data retrieved from Chembl)
match=as.vector(sapply(colnames(DRUG_DISTANCES),function(x){
  string= paste("\\b",x,"\\b",sep="")
  matches=grep(string,DrugTargets$MOLECULE_NAME,ignore.case=TRUE) 
}))
DTData_cmap=melt(match)
colnames(DTData_cmap)=c("row_orig_ChData","Cmap_drugs")

DTData_cmap=cbind(DTData_cmap,DrugTargets[DTData_cmap$row_orig_ChData,])
DTData_cmap=DTData_cmap[which(DTData_cmap$MECHANISM_OF_ACTION!=""),]
DTData_cmap=DTData_cmap[which(DTData_cmap$TARGET_NAME !=""),]
DTData_cmap=DTData_cmap[,-c(1,3)]
colnames(DTData_cmap)[1]="Drugs"
head(DTData_cmap)
DTData_cmap=unique(DTData_cmap)
write.csv(DTData_cmap,"../DATA/ProcessedDTargets.csv",row.names = FALSE)