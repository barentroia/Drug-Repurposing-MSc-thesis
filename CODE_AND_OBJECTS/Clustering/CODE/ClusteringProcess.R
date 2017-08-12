#Clustering Processing - This script has the necessary internal functions to generate the drug communities using each of the 5 algorithms and perform the enrichment analysis. To be used with EnrichmentAnalysis.R
library(apcluster)
library(cluster)
library(parallel)

#obtainCommTableALGORITHM is a function that generates a dataframe with drugs and communities cid. There is one for each algorithm, as the drug communities generation process is different.

#enrichedclustersALGORITHM performs the enrichment analysis for the clsutering results of that algorithm


##PAM
obtainCommTablePAM=function(i){
  pamres=pam(as.dist(DRUG_DISTANCES),i,diss=TRUE)
  pam_df=data.frame(cID=integer(1309),DRUGS=character(1309))
  pam_df$DRUGS=names(pamres$clustering)
  pam_df$cID=pamres$clustering
  return(pam_df)
}

enrichedclustersPAM=function(i,signific,level,filtered){
  pam_table=obtainCommTablePAM(i)
  if (filtered==TRUE){
    pam_table=filteredCommunities(pam_table, DRUG_DISTANCES)
  }
  return(enrichedclusters(pam_table,signific,level,filtered))
}
##MCL
obtainCommTableMCL=function(infl){
  infl=infl*10
  results.list<- strsplit(readLines(paste("../DATA/MCL_OUTPUT/out.MCL_SimTab","_noairq.txt.I",infl,sep="")),"\t")
  mclres=data.frame(cID=integer(1309),DRUGS=character(1309))
  mclres$DRUGS=unlist(results.list)
  j=1
  for (clusters in results.list){
    mclres$cID[mclres$DRUGS%in%clusters]=j
    j=j+1
  }
  return(mclres)
}

enrichedclustersMCL=function(infl,signific,level,filtered){
  mcl_table=obtainCommTableMCL(infl)
  if (filtered==TRUE){
    mcl_table=filteredCommunities(mcl_table, DRUG_DISTANCES)
  }
  res=enrichedclusters(mcl_table,signific,level,filtered)
  return(res)
}

##AP
#For comparison (note that I didnt remove any drugs.)
obtainCommTableAP=function(pref){
  apres=apcluster(1-DRUG_DISTANCES,details=TRUE,seed=10,p=pref) 
  Drugs=colnames(DRUG_DISTANCES)
  my_DRUG_COMMUNITIES=data.frame("cID"=numeric(length(Drugs)),"DRUGS"=Drugs)
  #Populating the matrix with the communities corresponding to each drug
  for (i in 1:length(apres)){
    my_DRUG_COMMUNITIES$cID[which(my_DRUG_COMMUNITIES$DRUGS%in%names(apres[[i]]))]=i
  }
  return(my_DRUG_COMMUNITIES)
}


enrichedclustersAP=function(p,signific,level,filtered){
  ap_table=obtainCommTableAP(p)
  if (filtered==TRUE){
    ap_table=filteredCommunities(ap_table, DRUG_DISTANCES)
  }
  return(enrichedclusters(ap_table,signific,level,filtered))
}

##HC
obtainCommTableHC=function(i){
  hclusters <- hclust(as.dist(DRUG_DISTANCES),method="complete")
  hcres <- cutree(hclusters, i)
  HC_df=data.frame(cID=integer(1309),DRUGS=character(1309))
  HC_df$DRUGS=names(hcres)
  HC_df$cID=hcres
  return(HC_df)
}

enrichedclustersHC=function(i,signific,level,filtered){
  hc_tab=obtainCommTableHC(i)
  if (filtered==TRUE){
    hc_tab=filteredCommunities(hc_tab, DRUG_DISTANCES)
  }
  return(enrichedclusters(hc_tab,signific,level,filtered))
}

##ClusterOne

obtainCommTableCL=function(dens){
  results.list<- strsplit(readLines(paste("../DATA/ClusterOne_Output/ClusterOne_d",dens,".txt",sep="")),"\t")
  nr=length(unlist(results.list))
  clres=data.frame(cID=integer(nr),DRUGS=character(nr))
  clres$DRUGS=unlist(results.list)
  repetitions=sapply(results.list,FUN=length)
  clres$cID=rep(c(1:length(results.list)),repetitions)
  
  return(clres)
}

enrichedclustersCL=function(dens,signific,level,filtered){
  cl_table=obtainCommTableCL(dens)
  if (filtered==TRUE){
    cl_table=filteredCommunities(cl_table, DRUG_DISTANCES)
  }
  return(enrichedclusters(cl_table,signific,level,filtered))
}

##Filtering done by Iorio

#' removeCol returns TRUE or FALSE depending on whether the drug column should be removed or not
#'
#' @param column vector with the distances between the drug under study and all other drugs
#' @return TRUE or FALSE
removecol=function(column){
  #note that in the drug column, the distance to itself is 0
  if (length(which(column<0.8))==1){
    return(TRUE) 
  }else { return(FALSE)}
}

filteredCommunities=function(tab,DRUG_DISTANCES){
  #DrugsToRemove will contain all the drugs that do not meet the threshold: they don't have a distance to other drugs in the same community smaller than 0.8
  DrugsToRemove=vector()
  for (i in 1:length(unique(tab$cID))){
    drugs=as.character(tab$DRUGS[which(tab$cID==i)])
    if (length(drugs)>1){
      Com_Matrix=DRUG_DISTANCES[which(rownames(DRUG_DISTANCES)%in%drugs),which(colnames(DRUG_DISTANCES)%in%drugs)]
      Drugs=apply(Com_Matrix,2,removecol)
      DrugsToRemove=c(DrugsToRemove,names(which(Drugs==TRUE)))
    }
    else{DrugsToRemove=c(DrugsToRemove,drugs)}
  }
  return(tab[which(!tab$DRUGS%in%DrugsToRemove),])
}

#This does the enrichment analysis for random clusters - the median values for the enrichedclusters etc are taken. 100 random clustering results are generated
randomParallel=function(tab,signific,level,filtered=FALSE){
cl=makeCluster(8,type = "FORK")
clusterExport(cl,varlist=ls(),envir=environment())
clusterEvalQ(cl, library(cluster))
clusterEvalQ(cl, library(mclust))
res_random=parSapply(cl,c(1:100),EnrichmentRandomClusters,tab=tab,signific=signific,level=level)
stopCluster(cl)
res=data.frame(K=quantile(as.numeric(res_random["K",]),0.50),EnrichedClusters=quantile(as.numeric(res_random["EnrichedClusters",]),0.50),Percentage=quantile(as.numeric(res_random["Percentage",]),0.50),EnrichedClustersDT=quantile(as.numeric(res_random["EnrichedClustersDT",]),0.50),PercentageDT=quantile(as.numeric(res_random["PercentageDT",]),0.50),Pval_ATC=quantile(as.numeric(res_random["Pval_ATC",]),0.50),Pval_DT=quantile(as.numeric(res_random["Pval_DT",]),0.50))
return(res)}
