#PAM
source("EnrichmentAnalysis.R")
source("ClusteringProcess.R")
#function to perform enrichment analysis of random communities
enrichedRandclustersPAM=function(i,signific,level,filtered){
  pam_table=obtainCommTablePAM(i)
  if (filtered==TRUE){
    pam_table=filteredCommunities(pam_table, DRUG_DISTANCES)
  }
  res_random=randomParallel(pam_table,signific,level)
  return(res_random)
}

ks=seq(2,400,by=4)
l=4
Results_PAM_koptimiz=sapply(ks,enrichedRandclustersPAM,signific=0.05,level=l,filtered=FALSE) 
Results_PAM_koptimiz=as.matrix(Results_PAM_koptimiz)
write.csv(Results_PAM_koptimiz,paste("../DATA/ParameterOptimization/Rand_PAM_koptimiz",l,".csv",sep=""))