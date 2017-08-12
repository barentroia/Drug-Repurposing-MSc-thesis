#AP Clustering - obtain 150 clustering results, by varying p parameter (preference)

library(apcluster)
source("EnrichmentAnalysis.R")
source("ClusteringProcess.R")

#various p
ps=c(-30,-10,-5,-2.2,-2,-1.5,-1,-0.75,-0.5,seq(-0.45,0.25,by=0.005))
lev=c(1,3,4,5) # several ATC levels

#non filtered - produce a table with information about each clustering result. 
for (l in lev){
Results_AP_koptimiz=sapply(ps,enrichedclustersAP,signific=0.05,level=l,filtered=FALSE) 
Results_AP_koptimiz=as.matrix(Results_AP_koptimiz)
write.csv(Results_AP_koptimiz,paste("../DATA/ParameterOptimization/Results_AP_koptimiz",l,".csv",sep=""))
}

#filtered - produce a table with information about each clustering result. 
for (l in lev){
  Results_AP_koptimiz=sapply(ps,enrichedclustersAP,signific=0.05,level=l,filtered=TRUE) 
  Results_AP_koptimiz=as.matrix(Results_AP_koptimiz)
  write.csv(Results_AP_koptimiz,paste("../DATA/ParameterOptimization/Results_AP_koptimiz_filtered",l,".csv",sep=""))
}
