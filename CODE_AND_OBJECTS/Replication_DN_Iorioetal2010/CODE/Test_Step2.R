#Test - confirm that PRLs produced are the same

print("Retrieving the distances matrix made available in Iorio 2015 paper")
load("../DATA/DRUG_DISTANCES.ro")
print("Retrieving the matrix with the computed distances")
my_Dist=read.table("../OUTPUT/DrugDistances.txt",check.names = FALSE,stringsAsFactors = FALSE,header=TRUE)
print("Reordering the columns of the distance matrix produced so that it can be compared with the distance matrix made available in Iorio 2015 paper")
#colnames(DRUG_DISTANCES)[which(!colnames(DRUG_DISTANCES)%in%colnames(my_Dist))] - yields 1
colnames(my_Dist)[which(!colnames(my_Dist)%in%colnames(DRUG_DISTANCES))]=colnames(DRUG_DISTANCES)[which(!colnames(DRUG_DISTANCES)%in%colnames(my_Dist))]
rownames(my_Dist)[which(!rownames(my_Dist)%in%rownames(DRUG_DISTANCES))]=rownames(DRUG_DISTANCES)[which(!rownames(DRUG_DISTANCES)%in%rownames(my_Dist))]

print(paste("Distance matrices have the same content (with 5e-05 tolerance):",all.equal(my_Dist,as.data.frame(DRUG_DISTANCES),tolerance=5e-05),sep=""))
