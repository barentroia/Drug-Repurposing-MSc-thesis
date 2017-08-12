source("EnrichmentAnalysis.R")
source("ClusteringProcess.R")

#AP
#function to perform enrichment analysis of random communities
enrichedRandclustersAP=function(p,signific,level,filtered){
  ap_table=obtainCommTableAP(p)
  if (filtered==TRUE){
    ap_table=filteredCommunities(ap_table, DRUG_DISTANCES)
  }
  res_random=randomParallel(ap_table,0.05,4)
  return(res_random)
}

ps=c(-30,-10,-5,-2.2,-2,-1.5,-1,-0.75,-0.5,seq(-0.45,0.25,by=0.005))
Results_AP_koptimiz=sapply(ps,enrichedRandclustersAP,signific=0.05,level=4,filtered=FALSE) 
Results_AP_koptimiz=as.matrix(Results_AP_koptimiz)
write.csv(Results_AP_koptimiz,paste("../DATA/ParameterOptimization/Rand_AP_koptimiz",4,".csv",sep=""))