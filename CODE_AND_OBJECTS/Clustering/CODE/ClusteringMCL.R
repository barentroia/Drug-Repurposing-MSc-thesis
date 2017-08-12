#MCL results -obtain 30 clustering results, by varying I parameter (inflation)

source("EnrichmentAnalysis.R")
source("ClusteringProcess.R")

#various i
infls=seq(1.4,7.0,0.2)
lev=c(1,3,4,5)

#non filtered - produce a table with information about each clustering result.
for (l in lev){
  Results_MCL_koptimiz=sapply(infls,enrichedclustersMCL,signific=0.05,level=l,filtered=FALSE) 
  Results_MCL_koptimiz=as.matrix(Results_MCL_koptimiz)
  write.csv(Results_MCL_koptimiz,paste("../DATA/ParameterOptimization/Results_MCL_koptimiz",l,".csv",sep=""))
}
#filtered - produce a table with information about each clustering result. 
for (l in lev){
  Results_MCL_koptimiz=sapply(infls,enrichedclustersMCL,signific=0.05,level=l,filtered=TRUE) 
  Results_MCL_koptimiz=as.matrix(Results_MCL_koptimiz)
  write.csv(Results_MCL_koptimiz,paste("../DATA/ParameterOptimization/Results_MCL_koptimiz_filtered",l,".csv",sep=""))
}

