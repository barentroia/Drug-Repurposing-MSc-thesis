#Replicating the drug network in Iorio 2010
## Step 1:  Build the PRLS
## Using the ConnectivityMap package 

packages <- c("ConnectivityMap")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library("ConnectivityMap")

print("Retrieving ConnectivityMap data. This takes 45 sec")
data(rankMatrix) # rankMatrix is a matrix with the gene expression profiles of each of the 6100 instances (i.e treatments of a certain cell line, with a certain drug, duration, concentration, described in instances matrix)
data(instances) # instances is a matrix with the description of the treatment in each instance(column in rankMatrix)

#Note: instances$cmap_name are used as drugs identifiers 
print("Retrieving the instances for each drug and saving them in GEPs of each drug folder. This takes approximately 3.5 minutes")
#Using the info in instances table, generate several subtables of rankMatrix, each one referring to all the instances belonging to the same drug 
#This for loop takes approximately 3 minutes. It genereates 1309 files, one for each drug and each one contains exactly one table. This table has as many columns as instances available in CMap for that drug
for (drug in unique(instances$cmap_name)){
  drug_instances=rownames(instances)[which(instances$cmap_name==drug)] # These are the identifiers of each instance belonging to the same drug
  
  #if statement so that a data frame is obtained and not a vector, regardless of how many instances each drug has
  if(length(drug_instances)>1){
  drug_GEPs=rankMatrix[,drug_instances]
  }else {
  drug_GEPs=data.frame(rankMatrix[,drug_instances],row.names=rownames(rankMatrix))
  }
  
  #some drugs would yield errors when saving the data tables, so I had to change slightly the names
  if(startsWith(drug,"\"")){
    drug = substring(drug, 2, nchar(drug)-1)
  }else if (grepl("/",drug)){
    drug=gsub("/", "", drug)
  }
  write.table(drug_GEPs, file=paste("GEPs of each drug/",drug,".out",sep=""), row.name=TRUE,col.names=TRUE)
}


#Functions to retrieve the PRLs for each drug : Spearmans Footrule, Borda Merging Function and Kruskal algorithm

#' SpeamansFootrule computes the distance between two ranked lists (columns i and j in rankMatrix)
#'
#' @param i The number of the column where one of the ranked Lists to be compared is
#' @param j The number of the column wherethe other ranked List is
#' @param rankMatrix the matrix to which the ranked Lists to be compared belong to. - retireved from ConnectivityMap package. Each column refers to one instance, each entry to a gene
#' @return dist Spearman's Footrule distance
SpearmansFootrule_ij<-function(i,j,rankMatrix){
  dist=sum(abs(rankMatrix[,i]-rankMatrix[,j]))
  return(dist)
}

#Vectorized the function so that I could use outer in Kruskal algorithm and avoid 2 for loops
SpearmansFootrule=Vectorize(SpearmansFootrule_ij,vectorize.args=list("i","j"))

  

#' Borda Merging Function - merges two ranked lists (GEPa and GEPb) GEP stands for gene expression profile
#'
#' @param GEPa One of the 2 gene expression profiles (ranked list) to be merged
#' @param GEPb The other gene expression profile (ranked list)
#' @return mergedGEP The result of merging the two GEPs
BordaMergingFunction=function(GEPa,GEPb){
  ps=GEPa+GEPb
  mergedGEP=rank(ps,ties.method="first") # In case of ties, the gene that is first (following an alphanumeric order) - Iorio used the same, as it yields the same results
  return(mergedGEP)
}



#' Kruskal Algorithm - to obtain a single ranked list from a set of ranked lists belonging to the same drug. This set is the input (GEPs)
#' Note that during the implementation, the initial set of lists is going to be modified, until only one ranked list is left
#'
#' @param GEPs The gene expression profiles that belong to the same drug and that are going to be merged
#' @return The result of merging the two GEPs
KruskalAlg=function(GEPs){
  n=length(GEPs)
  while (n>1){
    distmat=outer(1:n,1:n,SpearmansFootrule,rankMatrix=GEPs) # computes the distance between every pair of GEPs in the table provided, according to Spearmans Footrule
    #finds the two GEPs in the table which have the smallest distance
    min_ij=min(distmat[upper.tri(distmat)]) # note that the diagonal is always zero and distance is symmetrical
    i=which(distmat == min_ij, arr.ind = TRUE)[[1,1]]
    j=which(distmat == min_ij, arr.ind = TRUE)[[1,2]]
    #merges the most similar GEPs
    y=BordaMergingFunction(GEPs[,i],GEPs[,j])
    #replaces the most similar GEPs by the merged GEP
    GEPs[,i]=NULL
    GEPs[,j]=NULL
    GEPs=cbind(GEPs,y)
    n=length(GEPs)
  }
  
  #Some tidying before yielding the final PRL: PRL is a list with the genes ordered by their summarised differential expression
  GEPs$genes=rownames(GEPs)
  GEPs=GEPs[order(GEPs[,1]),]
  GEPs[,1]=NULL
  rownames(GEPs)=NULL
  return(GEPs)
}

#' retriveInstancesFile -  returns the data frame with the instances that refer to the drug provided
#'
#' @param drug The drug whose instances file we want to retrieve
#' @return drug_instances The data frame with the drug instances
retrieveInstancesFile=function(drug){
  #reminder: some drugs had names that caused problems when producing the files. So to find their files, the names have to be modified again
  if(startsWith(drug,"\"")){
    drug = substring(drug, 2, nchar(drug)-1)
  }else if (grepl("/",drug)){
    drug=gsub("/", "", drug)
  }
  drug_instances= read.table(paste("GEPs of each drug/",drug,".out",sep=""), header=TRUE)
  return(drug_instances)
}

#' generatePRLfromFile - generates the PRL (my_PRL) from the data file with the instances belonging to the same drug
#'
#' @param drug The drug whose PRL we want to generate
#' @return my_PRL The PRL that was generated
generatePRLfromFile=function(drug){
  drug_instances=retrieveInstancesFile(drug)
  my_PRL=KruskalAlg(drug_instances)
  
  #PRL should have the drug name as the column name to make it easier when building the PRL matrix
  if(startsWith(drug,"\"")){
    drug = substring(drug, 2, nchar(drug)-1)
  }
  drug=gsub(" ", "_", drug)
  colnames(my_PRL)=drug
  return(my_PRL)
}

#'buildPRLMatrix- builds the final matrix (PRLs), with all the PRLs built for each drug with generatePRLfromFile
#'
#' @param drugs The drugs that are represented in CMap drugs = unique(instances$cmap_name)
#' @return PRLs The DRUG_PRL matrix: each column represents a drug, each row a gene. Entries represent the ranks of the genes in the PRL of that drug
buildPRLMatrix=function(drugs){
  PRLs=data.frame(matrix(NA, nrow = 22283, ncol = 0))
  for (drug in drugs){
    #each PRL is a column with rownames= ranks, and column name = name of the drug (according to DRUG_PRLs standards)
    PRL=generatePRLfromFile(drug)
    #add PRL to PRLs data frame
    PRLs=cbind(PRLs,PRL)
  }
  PRLs=PRLs[ , order(colnames(PRLs))]
  return(PRLs)
}

print("Obtaining the PRLs of each drug. This takes approximately 9 minutes")
#Final step - Build my DRUG PRLs matrix
#Takes approximately 9 minutes
my_DRUG_PRLs=buildPRLMatrix(unique(instances$cmap_name))
print("Saving the matrix with all PRLs in Output folder. This takes 1 min")
write.table(my_DRUG_PRLs,"../OUTPUT/DrugPRLs.txt")