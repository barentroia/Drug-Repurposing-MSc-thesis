source("EnrichmentAnalysis.R")
source("ClusteringProcess.R")

#function to perform enrichment analysis of random communities
enrichedRandclustersCL=function(d,signific,level,filtered){
  cl_table=obtainCommTableCL(d)
  if (filtered==TRUE){
    cl_table=filteredCommunities(cl_table, DRUG_DISTANCES)
  }
  res_random=randomParallel(cl_table,signific,level)
  return(res_random)
}

dens = seq(0.0125,0.7,by=0.01)
l=4
Results_CL_koptimiz=sapply(dens,enrichedRandclustersCL,signific=0.05,level=l,filtered=FALSE) 
Results_CL_koptimiz=as.matrix(Results_CL_koptimiz)
write.csv(Results_CL_koptimiz,paste("../DATA/ParameterOptimization/Rand_CL_koptimiz",l,".csv",sep=""))