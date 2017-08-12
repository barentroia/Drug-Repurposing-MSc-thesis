source("EnrichmentAnalysis.R")
source("ClusteringProcess.R")

#function to perform enrichment analysis of random communities
enrichedRandclustersMCL=function(infl,signific,level,filtered){
  mcl_table=obtainCommTableMCL(infl)
  if (filtered==TRUE){
    mcl_table=filteredCommunities(mcl_table, DRUG_DISTANCES)
  }
  res=randomParallel(mcl_table,signific,level)
  return(res)
}

infls=seq(1.4,6.0,0.1)

Results_MCL_koptimiz=sapply(infls,enrichedRandclustersMCL,signific=0.05,level=4,filtered=FALSE) 
Results_MCL_koptimiz=as.matrix(Results_MCL_koptimiz)
write.csv(Results_MCL_koptimiz,paste("../DATA/ParameterOptimization/Rand_MCL_koptimiz",4,".csv",sep=""))