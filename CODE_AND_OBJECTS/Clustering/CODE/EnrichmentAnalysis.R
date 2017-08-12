#Enrichment Analysis 

#Reference data
MergedATCcodes=read.csv("../DATA/ProcessedATCcodes.csv",stringsAsFactors = FALSE,header=TRUE)
#MergedATCcodes$Codes=substring(MergedATCcodes$Codes,1,5)
#MergedATCcodes=unique(MergedATCcodes)

DrugTargets=read.csv("../DATA/ProcessedDTargets.csv",header=TRUE,stringsAsFactors = FALSE)
#DrugPairs=read.csv("../DATA/DrugPairsInfo.csv",stringsAsFactors = FALSE,header=TRUE)

#Distances between drugs
load("../../Replication_DN_Iorioetal2010/DATA/DRUG_DISTANCES.ro")
load("../../Replication_DN_Iorioetal2010/DATA/DRUG_COMMUNITIES.ro")
library(parallel)

Iorio_AP_DComm=read.table("../DATA/AP_nofilter.txt",header=TRUE) # Reference communities

#converting AP Iorio communities into vector
apres_vect=Iorio_AP_DComm$cID
names(apres_vect)=Iorio_AP_DComm$DRUGS

#obtaining the filtered communities
apres_vect_filtered=DRUG_COMMUNITIES$cID
names(apres_vect_filtered)=DRUG_COMMUNITIES$DRUGS
removedDrugs=names(apres_vect)[!names(apres_vect)%in%names(apres_vect_filtered)]
imaginaryComm=seq(4000,5400,by=1)
singleClusters=imaginaryComm[1:length(removedDrugs)]
names(singleClusters)=removedDrugs
apres_vect_filtered=c(apres_vect_filtered,singleClusters)

library(mclust)
library(cluster)

#Functions - Enrichment analysis for each community

#Fisher Test
#Input: atc - atc/MoA of interest
#       atc_codes_in_com - vector with the atcs/MoAs in the community
#       RefTable - either MergedATCCodes or DrugTargets
#       attribute - either "Codes" for ATCs, or "MECHANISM_OF_ACTION" for MoA
#Output: pvalue that assess the significance of the overrepresentation of that atc in the community
fishTest_in_com=function(atc,atc_codes_in_com,RefTable,attribute){
  
  atc_code_occur_in_com=length(atc_codes_in_com$Drugs[which(atc_codes_in_com[,attribute]==atc)])
  atc_code_occur_outside=length(RefTable$Drugs[which(RefTable[,attribute]==atc)])-atc_code_occur_in_com
  
  other_atc_code_in_com=length(unique(atc_codes_in_com$Drugs))-atc_code_occur_in_com #only drugs with known ATC are counted
  
  other_atc_code_outside = length(unique(RefTable$Drugs))-atc_code_occur_in_com- atc_code_occur_outside - other_atc_code_in_com
  
  cmatrix=matrix(c(atc_code_occur_in_com,other_atc_code_in_com,atc_code_occur_outside,other_atc_code_outside),ncol=2)
  
  test=fisher.test(cmatrix,alternative = "greater")
  return(test$p.value)
}
#Retrieves the information needed to be used in fishTest_in_com function and uses it to perform the enrichement tests for the atc/moa in a community
#Input: i - cID of the community (community identifier)
#       DRUG_COMMUNITIES : data frame with two columns, one for the cIDs and the other for the corresponding drugs
#       RefTable - either MergedATCCodes or DrugTargets
#       attribute - either "Codes" for ATCs, or "MECHANISM_OF_ACTION" for MoA
#Output:p values of the enrichment tests for each atc of the community
performEnrichTestfinCommunity=function(i,DRUG_COMMUNITIES,RefTable,attribute){
  #retrieve the drugs in the community
  drugs_in_com=DRUG_COMMUNITIES$DRUGS[which(DRUG_COMMUNITIES$cID==i)] #1
  
  #retrieve their ATC codes
  atc_codes_in_com=RefTable[RefTable$Drugs%in%drugs_in_com,]
  
  #enrichment test is only done for atc codes that correspond to at least two drugs
  
  frequency_codes=table(atc_codes_in_com[,attribute])
  freq_codes=names(frequency_codes)[which(frequency_codes>1)]
  
  
  if(length(freq_codes)!=0){
    res=sapply(freq_codes,function(atc){fishTest_in_com(atc,atc_codes_in_com,RefTable,attribute)})
  }else{
    res=NA 
    names(res)=NA
  }
  return(res)
}

#Performs the enrichment analysis of all the communties in a clustering result, using the above functions
#Input: DRUG_COMMUNITIES : data frame with two columns, one for the cIDs and the other for the corresponding drugs
#       RefTable - either MergedATCCodes or DrugTargets
#       attribute - either "Codes" for ATCs, or "MECHANISM_OF_ACTION" for MoA
#Output:dataframe with the atc/moa in each community and the result of the fisher test

enrichmentClustering=function(DRUG_COMMUNITIES,RefTable,attribute){
  communities=unique(DRUG_COMMUNITIES$cID)
  list_res= lapply(communities,function(x){performEnrichTestfinCommunity(x, DRUG_COMMUNITIES,RefTable,attribute)})
  if (is.list(list_res)){
    ResultsEnrichment=data.frame(cID=rep(communities,sapply(list_res, FUN=length)),attribute=as.vector(names(unlist(list_res))),Pvalues=as.numeric(unlist(list_res)))
  }else{
    ResultsEnrichment=data.frame(cID=rep(communities,sapply(list_res, FUN=length)),attribute=as.vector(names(list_res)),Pvalues=as.numeric(list_res))
  }
  
  
  colnames(ResultsEnrichment)=c("cID",attribute,"Pvalues")
  
  
  #Correction for multiple testing
  ResultsEnrichment$Pvalues=p.adjust(ResultsEnrichment$Pvalues, method = "BH")
  return(ResultsEnrichment)
}

# Returns a data frame with information about the clustering result, namely enrichment analysis of its communities in ATC and MoA
#Input: tab - data frame with two columns, one for the cIDs and the other for the corresponding drugs
#       signif - the significance level for the fisher tests
#       level  - the ATC level that should be considered :  we use 4 -> to correspond to the 3rd level
#       filtered - TRUE or FALSE depending on whether we want to test the communities after they have been filtered or not
#Output:Data frame with : number of clusters (K),Drugs that were not filtered out, No. Enriched clusters in ATC/MoA, Percentage of Enr ATC/MoA clusters out of total clusters, Fisher test pvalue on whether drugs with similar MoA/DT tend to be in the same communities,ARI compared to the reference communtiies,Overlap with the enriched communities found)

enrichedclusters=function(tab,signific,level,filtered){
  MergedATCcodes$Codes=substring(MergedATCcodes$Codes,1,level)
  MergedATCcodes=unique(MergedATCcodes)
  
  #Performing the enrichment of Iorio's communities to then compared with the alternative clustering results 
  en_ATC_Iorio=enrichmentClustering(Iorio_AP_DComm,MergedATCcodes,"Codes")
  en_ATC_Iorio=as.character(unique(en_ATC_Iorio$Codes[which(en_ATC_Iorio$Pvalues<0.05)]))
  
  en_DT_Iorio=enrichmentClustering(Iorio_AP_DComm,DrugTargets,"MECHANISM_OF_ACTION")
  en_DT_Iorio=as.character(unique(en_DT_Iorio$MECHANISM_OF_ACTION[which(en_DT_Iorio$Pvalues<0.05)]))
  
  #Enrichment analysis of the clustering results
  en=enrichmentClustering(tab,MergedATCcodes,"Codes")
  enr_clusters=length(unique(en$cID[which(en$Pvalues<signific)]))
  
  enr_ATC=as.character(unique(en$Codes[which(en$Pvalues<0.05)]))
  intersect_Iorio=intersect(en_ATC_Iorio,enr_ATC)
  percentage_overlap_ATC=length(intersect_Iorio)/length(en_ATC_Iorio)
  
  en_DT=enrichmentClustering(tab,DrugTargets,"MECHANISM_OF_ACTION")
  enr_clusters_DT=length(unique(en_DT$cID[which(en_DT$Pvalues<signific)]))
  
  enr_DT=as.character(unique(en_DT$MECHANISM_OF_ACTION[which(en_DT$Pvalues<0.05)]))
  intersect_DT_Iorio=intersect(en_DT_Iorio,enr_DT)
  percentage_overlap_DT=length(intersect_DT_Iorio)/length(en_DT_Iorio)
  
  #Compute other info about the clustering
  RemovedDrugs=1309-length(unique(tab$DRUGS))
  IncludedDrugs=1309-RemovedDrugs-length(which(table(tab$cID)==1))
  totalClusters=length(which(table(tab$cID)>1))
  
  percentage=round(enr_clusters*100/totalClusters,2)
  percentage_DT=round(enr_clusters_DT*100/totalClusters,2)
  
  #1 Fisher Test
  DrugPairs=read.csv(paste("ATC CODES/DrugPairsInfolevel",level,".csv",sep=""),stringsAsFactors = FALSE,header=TRUE)
  DrugPairs_k=addCommunitiestoDrugPairs(DrugPairs,tab)
  pvalATC=runFisherTestDrugPairs(DrugPairs_k,"SameATC")
  pvalDT=runFisherTestDrugPairs(DrugPairs_k,"SameDrugTarget")
  vect_res=as.numeric(tab$cID)
  names(vect_res)=tab$DRUGS
  
  
  #Adjusted Rand Index
  if (filtered==TRUE){
    removedDrugs=names(apres_vect)[!names(apres_vect)%in%names(vect_res)]
    singleClusters=imaginaryComm[1:length(removedDrugs)]
    names(singleClusters)=removedDrugs
    vect_res=c(vect_res,singleClusters)
    vect_res=vect_res[match(names(apres_vect_filtered),names(vect_res))]
    adj=adjustedRandIndex(apres_vect_filtered,vect_res)
  }else{
    vect_res=vect_res[match(names(apres_vect),names(vect_res))]
    adj=adjustedRandIndex(apres_vect,vect_res)
  }
  

  
  res=data.frame(K=length(unique(tab$cID)),"IncDrugs"= IncludedDrugs, "EnrichedClusters"=enr_clusters, "Percentage"=percentage, "EnrichedClustersDT"=enr_clusters_DT, "PercentageDT"=percentage_DT, "Pval_ATC"=pvalATC,"Pval_DT"=pvalDT,"ARIndex"=adj,"Overlap_ATC"=percentage_overlap_ATC,"Overlap_DT"=percentage_overlap_DT)
  return(res)
}

#Random Clustering
EnrichmentRandomClusters=function(x,tab,signific,level){
  my_DRUG_COMMUNITIES_random=tab
  my_DRUG_COMMUNITIES_random$cID=sample(my_DRUG_COMMUNITIES_random$cID)
  return(enrichedclusters(my_DRUG_COMMUNITIES_random,signific,level,FALSE))}

#Function - Only one global Fisher Exact Test: test whether drugs with similar MoA/DT tend to be in the same communities

#building a dataframe that for each pair of drugs asserts if they 1) share the same MoA/ATC 2) are in the same community
buildDrugPairsTable=function(DRUG_DISTANCES,MergedATCcodes,DrugTargets){
  #retrieve all drug pairs
  DrugPairs=DRUG_DISTANCES
  DrugPairs[upper.tri(DrugPairs)]=NA
  diag(DrugPairs)=NA
  DrugPairs=melt(DrugPairs, na.rm = TRUE)
  DrugPairs$value=NULL
  colnames(DrugPairs)=c("Drug_1","Drug_2")
  
  #Identify which drugs have the same ATC
  DrugPairs$SameATC=sapply(c(1:nrow(DrugPairs)),function(x){
    same=intersect(MergedATCcodes$Codes[which(MergedATCcodes$Drugs==DrugPairs$Drug_1[x])],MergedATCcodes$Codes[which(MergedATCcodes$Drugs==DrugPairs$Drug_2[x])])
    if (length(same)==0){return(FALSE)} else{return(TRUE)}})
  
  DrugsATC=unique(MergedATCcodes$Drugs)
  DrugPairs$SameATC[!DrugPairs$Drug_1%in%DrugsATC]="Not known"
  DrugPairs$SameATC[!DrugPairs$Drug_2%in%DrugsATC]="Not known"
  
  #Identify which drugs have the same MoA
  DrugPairs$SameDrugTarget=sapply(c(1:nrow(DrugPairs)),function(x){
    same=intersect(DrugTargets$MECHANISM_OF_ACTION[which(as.character(DrugTargets$Drugs)==DrugPairs$Drug_1[x])],DrugTargets$MECHANISM_OF_ACTION[which(as.character(DrugTargets$Drugs)==DrugPairs$Drug_2[x])])
    if (length(same)==0){return(FALSE)} else{return(TRUE)}})
  
  DrugsDT=unique(as.character(DrugTargets$Drugs))
  DrugPairs$SameDrugTarget[!DrugPairs$Drug_1%in%DrugsDT]="Not known"
  DrugPairs$SameDrugTarget[!DrugPairs$Drug_2%in%DrugsDT]="Not known"
  return(DrugPairs)}


#Attaches a column to the dataframe produced by buildDrugPairsTable indicating whether the drug pair is in the same community or not
addCommunitiestoDrugPairs=function(DrugPairs,my_DRUG_COMMUNITIES){
#Label each pair as belonging to the Same Community or Different

DrugPairs$comm1=my_DRUG_COMMUNITIES$cID[match(DrugPairs$Drug_1,my_DRUG_COMMUNITIES$DRUGS)]
DrugPairs$comm2=my_DRUG_COMMUNITIES$cID[match(DrugPairs$Drug_2,my_DRUG_COMMUNITIES$DRUGS)]
DrugPairs$SameComm=FALSE
DrugPairs$SameComm[which(DrugPairs$comm1==DrugPairs$comm2)]=TRUE
return(DrugPairs)}


#based on the DrugPairs table, built with the two functions above, performs a Fisher exact test to test whether drugs with similar MoA/DT tend to be in the same communities

runFisherTestDrugPairs=function(DrugPairs,attribute){
  sameComsameAtt=length(which(DrugPairs$SameComm==TRUE & DrugPairs[,attribute]==TRUE))
  sameComdifAtt=length(which(DrugPairs$SameComm==TRUE & DrugPairs[,attribute]==FALSE))
  difComsameAtt=length(which(DrugPairs$SameComm==FALSE & DrugPairs[,attribute]==TRUE))
  difComdifAtt=length(which(DrugPairs$SameComm==FALSE & DrugPairs[,attribute]==FALSE))
  
  cmatrix=matrix(c(sameComsameAtt,sameComdifAtt,difComsameAtt,difComdifAtt),ncol=2)
  test=fisher.test(cmatrix,alternative="greater")
  return(test$p.value)
}

