#PAM - Clustering obtain 20 clustering results, by varying k parameter (no clusters)
library(cluster)
source("EnrichmentAnalysis.R")
source("ClusteringProcess.R")

#various ks
ks=seq(2,400,by=2)
lev=c(1,3,5) #several ATC levels

#non filtered - produce a table with information about each clustering result.
for (l in lev){
  Results_PAM_koptimiz=sapply(ks,enrichedclustersPAM,signific=0.05,level=l,filtered=FALSE) 
  Results_PAM_koptimiz=as.matrix(Results_PAM_koptimiz)
  write.csv(Results_PAM_koptimiz,paste("../DATA/ParameterOptimization/Results_PAM_koptimiz",l,".csv",sep=""))
}

#filtered - produce a table with information about each clustering result. 

for (l in lev){
  Results_PAM_koptimiz=sapply(ks,enrichedclustersPAM,signific=0.05,level=l,filtered=TRUE) 
  Results_PAM_koptimiz=as.matrix(Results_PAM_koptimiz)
  write.csv(Results_PAM_koptimiz,paste("../DATA/ParameterOptimization/Results_PAM_koptimiz_filtered",l,".csv",sep=""))
}