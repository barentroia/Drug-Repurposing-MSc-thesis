#ClusterOne - Clustering - obtain 74 clustering results, by varying dens parameter

#code for console
#for dens in $(seq 0.0125 0.01 1); do java -jar cluster_one-1.0.jar CL_SimTab_betw0_1_filtered.txt -s 1 -d $dens > #"ClusterOne_Output/ClusterOne_d${dens}.txt"; done
source("EnrichmentAnalysis.R")
source("ClusteringProcess.R")

#various densities
dens = seq(0.0125,0.7425,by=0.01)

lev=c(1,3,4,5)

#non filtered - produce a table with information about each clustering result. 
for (l in lev){
  Results_CL_koptimiz=sapply(dens,enrichedclustersCL,signific=0.05,level=l,filtered=FALSE) 
  Results_CL_koptimiz=as.matrix(Results_CL_koptimiz)
  write.csv(Results_CL_koptimiz,paste("../DATA/ParameterOptimization/Results_CL_koptimiz",l,".csv",sep=""))
}

#filtered - produce a table with information about each clustering result. 
for (l in lev){
  Results_CL_koptimiz=sapply(dens,enrichedclustersCL,signific=0.05,level=l,filtered=TRUE) 
  Results_CL_koptimiz=as.matrix(Results_CL_koptimiz)
  write.csv(Results_CL_koptimiz,paste("../DATA/ParameterOptimization/Results_CL_koptimiz_filtered",l,".csv",sep=""))
}
