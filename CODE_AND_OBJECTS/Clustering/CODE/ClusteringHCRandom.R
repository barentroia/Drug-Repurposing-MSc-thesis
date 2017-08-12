source("EnrichmentAnalysis.R")
source("ClusteringProcess.R")

#function to perform enrichment analysis of random communities
enrichedRandclustersHC=function(i,signific,level,filtered){
  hc_tab=obtainCommTableHC(i)
  if (filtered==TRUE){
    hc_tab=filteredCommunities(hc_tab, DRUG_DISTANCES)
  }
  res_random=randomParallel(hc_tab,signific,level)
  return(res_random)
}

ks=seq(2,400,by=4)

Results_HC_koptimiz=sapply(ks,enrichedRandclustersHC,signific=0.05,level=4,filtered=FALSE)
Results_HC_koptimiz=as.matrix(Results_HC_koptimiz)
write.csv(Results_HC_koptimiz,paste("../DATA/ParameterOptimization/Rand_HC_koptimiz",4,".csv",sep=""))