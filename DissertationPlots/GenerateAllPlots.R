#Script to generate all of the plots in the dissertation

#Calling from DissertationPlots
cur.direct=getwd()
#Exploratory CMap Data Analysis
source("CMAP_expldataanalysisFP.R")
#Clustering Plots
source("../CODE_AND_OBJECTS/Clustering/CODE/FinalPlots.R")
setwd(cur.direct)
#Proximity of clusters
source("../CODE_AND_OBJECTS/QueryingDvD/CODE/InvestigatingClustersProximity.R")
setwd(cur.direct)