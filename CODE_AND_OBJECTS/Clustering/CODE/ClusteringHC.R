#HC - Clustering obtain 200 clustering results, by varying k parameter (no clusters)

source("EnrichmentAnalysis.R")
source("ClusteringProcess.R")

#various ks
ks=seq(2,400,by=2)
lev=c(1,3,5)

#non filtered - produce a table with information about each clustering result. 
for (l in lev){
  Results_HC_koptimiz=sapply(ks,enrichedclustersHC,signific=0.05,level=l,filtered=FALSE)
  Results_HC_koptimiz=as.matrix(Results_HC_koptimiz)
  write.csv(Results_HC_koptimiz,paste("../DATA/ParameterOptimization/testResults_HC_koptimiz",l,".csv",sep=""))
}

#filtered - produce a table with information about each clustering result. 
for (l in lev){
  Results_HC_koptimiz=sapply(ks,enrichedclustersHC,signific=0.05,level=l,filtered=TRUE)
  Results_HC_koptimiz=as.matrix(Results_HC_koptimiz)
  write.csv(Results_HC_koptimiz,paste("../DATA/ParameterOptimization/Results_HC_koptimiz_filtered",l,".csv",sep=""))
}