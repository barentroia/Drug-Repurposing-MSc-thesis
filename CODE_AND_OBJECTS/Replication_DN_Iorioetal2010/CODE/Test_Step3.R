#Test - confirm that drug communities found are the same

print("Retrieving the Communities dataframe made available in Iorio 2015 paper")
load("../DATA/DRUG_COMMUNITIES.ro")
print("Retrieving the dataframe with the communities found")
my_DRUG_COMMUNITIES=read.table("../OUTPUT/DrugCommunities.txt",header=TRUE,check.names = FALSE,stringsAsFactors = FALSE)
print("Reordering the columns of the communities dataframe produced so that it can be compared with the communities dataframe made available in Iorio 2015 paper")
#changing some ordering details so that it looks exactly the same as DRUG_COMMUNITIES
my_DRUG_COMMUNITIES=my_DRUG_COMMUNITIES[order(my_DRUG_COMMUNITIES$cID,my_DRUG_COMMUNITIES$DRUGS),]
my_DRUG_COMMUNITIES$cID=as.numeric(my_DRUG_COMMUNITIES$cID)
#Check obtained matrix is the same as DRUG_COMMUNITIES
print(paste("Communities dataframes have the same content:",identical(my_DRUG_COMMUNITIES,DRUG_COMMUNITIES),sep=""))

#TEST - Returns the number of communities that are different from DRUG_COMMUNITIES matrix provided by Iorio et al 2015 (if needed)
equalCommunities=function(apres){
  exemplars = names(apres@exemplars)
  nr=1
  eq=vector()
  DRUGS_1233=DRUG_COMMUNITIES$DRUGS # Only 1233 drugs are represented in the drug network
  for (i in exemplars){
    NR_COM=DRUG_COMMUNITIES$cID[which(DRUG_COMMUNITIES$DRUGS==i)] # cID of the examplar in Iorio's communities
    DRUG_COM=DRUG_COMMUNITIES$DRUGS[which(DRUG_COMMUNITIES$cID==NR_COM)] # extract the drugs in that community
    DRUG_AP=names(apres[[nr]]) # drugs in my equivalent community (obtained with my code)
    DRUG_AP=intersect(DRUG_AP,DRUGS_1233) # Since Iorio has left out some drugs, the test only checks if the common drugs are in the same community
    eq[nr]=setequal(DRUG_COM,DRUG_AP)
    nr=nr+1
  }
  print(which(eq==FALSE))
  print(paste("There were:", length(which(eq==FALSE)),"different communities."))
}