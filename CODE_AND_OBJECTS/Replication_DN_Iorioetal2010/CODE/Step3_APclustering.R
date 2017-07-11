#Replicating the drug network in Iorio 2010
## Step 3: Finding drug communities, using Affinity propagation clustering algorithm

##Clustering: from DRUG_DISTANCES to DRUG_COMMUNITIES

packages <- c("apcluster")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

#Obtain the final network
library(apcluster)
load("../DATA/DRUG_DISTANCES.ro")

print("Performing the clustering of the drugs, using Affinity Propagation algorithm and the drug distances computed previously. This takes 1 min.")
##This command outputs an APResult object with the found clusters with their exemplars and drugs in each comunity
apres=apcluster(1-DRUG_DISTANCES,details=TRUE,seed=10)


#After obtaining the communities, we remove drugs whose distances to other drugs in the same community are all larger than 0.8

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

#' buildCommunitiesMatrix rearranges the output of apcluster (apres) to a dataframe, so that it can be more easily manipulable and compared with the communities matrix made available by Iorio et al 2015
#'
#' @param apres results of applying apcluster to the drug distances matrix
#' @param DRUG_DISTANCES matrix with the distances between drugs
#' @return communities dataframe
buildCommunitiesMatrix=function(apres,DRUG_DISTANCES){
#DrugsToRemove will contain all the drugs that do not meet the threshold: they don't have a distance to other drugs in the same community smaller than 0.8
DrugsToRemove=vector()
for (i in 1:length(apres)){
  Com_Matrix=DRUG_DISTANCES[which(rownames(DRUG_DISTANCES)%in%names(apres[[i]])),which(colnames(DRUG_DISTANCES)%in%names(apres[[i]]))]
  Drugs=apply(Com_Matrix,2,removecol)
  DrugsToRemove=c(DrugsToRemove,names(which(Drugs==TRUE)))
}
#So.. we can obtain Drug_COMMUNITIES replicate by
Drugs=as.character(colnames(DRUG_DISTANCES)[which(!colnames(DRUG_DISTANCES)%in%DrugsToRemove)]) # drugs that meet the threshold
#DCOM_rep is my DRUG_COMMUNITIES matrix, obtained with my code
my_DRUG_COMMUNITIES=data.frame("cID"=numeric(length(Drugs)),"DRUGS"=Drugs)
#Populating the matrix with the communities corresponding to each drug
for (i in 1:length(apres)){
  my_DRUG_COMMUNITIES$cID[which(my_DRUG_COMMUNITIES$DRUGS%in%names(apres[[i]]))]=i
}
rownames(my_DRUG_COMMUNITIES)=my_DRUG_COMMUNITIES$DRUGS
my_DRUG_COMMUNITIES$DRUGS=as.character(my_DRUG_COMMUNITIES$DRUGS)
return(my_DRUG_COMMUNITIES)}

print("Rearranging the output of the clsutering to a dataframe that is more easily manipulable")
my_DRUG_COMMUNITIES=buildCommunitiesMatrix(apres,DRUG_DISTANCES)

print("Saving the drug communities matrix in Output folder")
#Final Output - my_DRUG_COMMUNITIES
write.table(my_DRUG_COMMUNITIES,"../OUTPUT/DrugCommunities.txt",sep="\t")

#Rich Clubs 1
examplars1_distance=DRUG_DISTANCES[names(apres@exemplars),names(apres@exemplars)]
apres_richClubs1=apcluster(1-examplars1_distance,details=TRUE)
my_richClubs1=buildCommunitiesMatrix(apres_richClubs1,examplars1_distance) # yields 10 communities, 9 of them are exactly the same as the ones on Iorio excel table,but one does not exist (community 2).

examplars2_distance=DRUG_DISTANCES[names(apres_richClubs1@exemplars),names(apres_richClubs1@exemplars)]
apres_richClubs2=apcluster(1-examplars2_distance,details=TRUE)
my_richClubs2=buildCommunitiesMatrix(apres_richClubs2,examplars2_distance) # similar, apart from the examplar of the tenth community that does not appear in first rich clubs

