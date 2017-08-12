#Final Plots
setwd("../CODE_AND_OBJECTS/Clustering/CODE/")
source("ClusteringProcess.R")
source("EnrichmentAnalysis.R")

#Loading needed packages
packages <- c("ggsci", "ggplot2","pheatmap","ggthemes","reshape2")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library("ggsci")
library("ggplot2")
library("gridExtra")
library("grid")
library("ggthemes")
library("pheatmap")
library("reshape2")

#standard theme for all plots

theme_Publication <- function(base_size=14, base_family="helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1.5)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#manual colorfill so that the algorithms are always assigned the same colors
myColors = c("royalblue1","gold","gray65","red3","olivedrab2")
names(myColors) = levels(as.factor(c("AP","HC","MCL","PAM","CL1")))
colScale = scale_colour_manual(values = myColors)
colfillScale=scale_fill_manual(values=myColors)

#Create plot function
compareEnrichedClusters=function(AlgComparis,compar,ylab){
  subsetPerc_ATC=AlgComparis[,c("K",compar,"algorithm")]
  subsetPerc_ATC=melt(subsetPerc_ATC,id=c("algorithm","K"))
  ggplot(data=subsetPerc_ATC, aes(x=K, y=value, colour=algorithm))+ 
    geom_line(size=1.08)+ geom_vline(xintercept=106,linetype="dashed")+labs(x="Number of Clusters",y = ylab)+labs(colour = "")+theme_bw()+xlim(0,400)+colScale+theme_Publication()}

#Create and save Plot
generatePlot=function(AlgComparis,attribute,ylab,filename){
  compareEnrichedClusters(AlgComparis,attribute,ylab)
  ggsave(paste("../../../DissertationPlots/Comparison_",attribute,filename,".png",sep=""),width=13,height=12,units = "cm")
}


##Algorithms - Retrieving the info of several clustering results and tidying it up into a dataframe

tidyParamOptimResTable=function(filename,alg){
  tab=read.csv(filename,header=TRUE)
  names=as.character(tab$X)
  tab=t(tab[,-1])
  colnames(tab)=names
  tab=as.data.frame(tab)
  tab$algorithm=alg
  return(tab)
}
#HC
HCsevKs_level1=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_HC_koptimiz1.csv","HC")
HCsevKs_level3=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_HC_koptimiz3.csv","HC")
HCsevKs_level4=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_HC_koptimiz4.csv","HC")
HCsevKs_level5=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_HC_koptimiz5.csv","HC")

listHC=list(Level1=HCsevKs_level1,Level2=HCsevKs_level3,Level3=HCsevKs_level4,Level4=HCsevKs_level5)

HCsevKs_level1_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_HC_koptimiz_filtered1.csv","HC")
HCsevKs_level3_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_HC_koptimiz_filtered3.csv","HC")
HCsevKs_level4_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_HC_koptimiz_filtered4.csv","HC")
HCsevKs_level5_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_HC_koptimiz_filtered5.csv","HC")

listHC_f=list(Level1=HCsevKs_level1_f,Level2=HCsevKs_level3_f,Level3=HCsevKs_level4_f,Level4=HCsevKs_level5_f)
#AP
APsevPs_level1=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_AP_koptimiz1.csv","AP")
APsevPs_level3=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_AP_koptimiz3.csv","AP")
APsevPs_level4=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_AP_koptimiz4.csv","AP")
APsevPs_level5=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_AP_koptimiz5.csv","AP")

listAP=list(Level1=APsevPs_level1,Level2=APsevPs_level3,Level3=APsevPs_level4,Level4=APsevPs_level5)

APsevPs_level1_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_AP_koptimiz_filtered1.csv","AP")
APsevPs_level3_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_AP_koptimiz_filtered3.csv","AP")
APsevPs_level4_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_AP_koptimiz_filtered4.csv","AP")
APsevPs_level5_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_AP_koptimiz_filtered5.csv","AP")
listAP_f=list(Level1=APsevPs_level1_f,Level2=APsevPs_level3_f,Level3=APsevPs_level4_f,Level4=APsevPs_level5_f)

#MCL
MCLsevIs_level1=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_MCL_koptimiz1.csv","MCL")
MCLsevIs_level3=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_MCL_koptimiz3.csv","MCL")
MCLsevIs_level4=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_MCL_koptimiz4.csv","MCL")
MCLsevIs_level5=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_MCL_koptimiz5.csv","MCL")

listMCL=list(Level1=MCLsevIs_level1,Level2=MCLsevIs_level3,Level3=MCLsevIs_level4,Level4=MCLsevIs_level5)

MCLsevIs_level1_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_MCL_koptimiz_filtered1.csv","MCL")
MCLsevIs_level3_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_MCL_koptimiz_filtered3.csv","MCL")
MCLsevIs_level4_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_MCL_koptimiz_filtered4.csv","MCL")
MCLsevIs_level5_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_MCL_koptimiz_filtered5.csv","MCL")

listMCL_f=list(Level1=MCLsevIs_level1_f,Level2=MCLsevIs_level3_f,Level3=MCLsevIs_level4_f,Level4=MCLsevIs_level5_f)

#PAM

PAMsevKs_level1=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_PAM_koptimiz1.csv","PAM")
PAMsevKs_level3=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_PAM_koptimiz3.csv","PAM")
PAMsevKs_level4=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_PAM_koptimiz4.csv","PAM")
PAMsevKs_level5=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_PAM_koptimiz5.csv","PAM")

listPAM=list(Level1=PAMsevKs_level1,Level2=PAMsevKs_level3,Level3=PAMsevKs_level4,Level4=PAMsevKs_level5)

PAMsevKs_level1_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_PAM_koptimiz_filtered1.csv","PAM")
PAMsevKs_level3_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_PAM_koptimiz_filtered3.csv","PAM")
PAMsevKs_level4_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_PAM_koptimiz_filtered4.csv","PAM")
PAMsevKs_level5_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_PAM_koptimiz_filtered5.csv","PAM")
listPAM_f=list(Level1=PAMsevKs_level1_f,Level2=PAMsevKs_level3_f,Level3=PAMsevKs_level4_f,Level4=PAMsevKs_level5_f)


#CLOne
CLsevDs_level1=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_CL_koptimiz1.csv","CL1")
CLsevDs_level3=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_CL_koptimiz3.csv","CL1")
CLsevDs_level4=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_CL_koptimiz4.csv","CL1")
CLsevDs_level5=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_CL_koptimiz5.csv","CL1")

CLsevDs_level4=CLsevDs_level4[which(CLsevDs_level4$IncDrugs>500),]
CLsevDs_level4_f=tidyParamOptimResTable("../DATA/ParameterOptimization/Results_CL_koptimiz_filtered4.csv","CL1")
CLsevDs_level4_f=CLsevDs_level4_f[which(CLsevDs_level4_f$IncDrugs>500),]

######Plots######

#Adjusted RandIndex - before filtering
Algorithms=rbind(HCsevKs_level4,APsevPs_level4,MCLsevIs_level4,PAMsevKs_level4,CLsevDs_level4)
compareEnrichedClusters(Algorithms,"ARIndex","Adjusted Rand Index")
generatePlot(Algorithms,"ARIndex","Adjusted Rand Index","ARIndexALLNFiltered")

#ATC and Drug Targets
#ATC and DT trends
generateMeltedATCDT=function(lev){
ATC_DTcomparis=data.frame()

algcomp=list(listHC[[lev]],listAP[[lev]],listPAM[[lev]])

for (df in algcomp){
  subset=df[,c("K","EnrichedClusters","EnrichedClustersDT","algorithm")]
  ATC_DTcomparis=rbind(ATC_DTcomparis,subset)
}

ks=intersect(listHC[[lev]]$K,listPAM[[lev]]$K)
ks=c(ks,intersect(ks,listAP[[lev]]$K))
medianEnr=function(k){
  enrATC=median(c(ATC_DTcomparis$EnrichedClusters[which(ATC_DTcomparis$K==k)],1))
  enrDT=median(c(ATC_DTcomparis$EnrichedClustersDT[which(ATC_DTcomparis$K==k)],1))
  return(list(ATC=enrATC,DT=enrDT))}
values=as.data.frame(t(as.data.frame(sapply(ks,medianEnr))))
medians=data.frame(K=ks,"ATC codes (3rd level)"=as.numeric(values$ATC),"Mode of Action"=as.numeric(values$DT))
subset_m=melt(medians,id="K")
return(subset_m)
}
subset_1=generateMeltedATCDT(1)
subset_2=generateMeltedATCDT(2)
subset_3=generateMeltedATCDT(3)
subset_4=generateMeltedATCDT(4)

ggplot(data=subset_3, aes(x=K, y=value, colour=variable)) +
  geom_line(size=1.08)+labs(x="Number of Clusters",y = "Number of Enriched Communities")+xlim(0,400)+labs(colour = "")+theme_bw()+scale_color_jco(labels = c("ATC Code (3rd level)", "Mode of Action"))+theme_Publication()
ggsave("../../../DissertationPlots/Comparison_ATC3rd_DT.png",width=13,height=12,units = "cm")


#Heatmap - clustering results at 106
MergedATCcodes=read.csv("../DATA/ProcessedATCcodes.csv",stringsAsFactors = FALSE,header=TRUE)
DrugTargets=read.csv("../DATA/ProcessedDTargets.csv",header=TRUE,stringsAsFactors = FALSE)

obtainEnrichedSet=function(tab,attribute,RefTable,level=5,signific=0.05){
  if (attribute =="Codes"){
    RefTable$Codes=substring(RefTable$Codes,1,level)
    RefTable=unique(RefTable)
  }
  res=enrichmentClustering(tab,RefTable,attribute)
  enriched=res[which(res$Pvalues<signific),][,attribute]
  enriched=as.character(unique(enriched))
  return(enriched)
}

AP=obtainCommTableAP(pref=0.0681)
#AP_f=filteredCommunities(AP,DRUG_DISTANCES)

MCL=obtainCommTableMCL(2.2)
#MCL_f=filteredCommunities(MCL,DRUG_DISTANCES)

PAM=obtainCommTablePAM(106)
#PAM_f=filteredCommunities(PAM,DRUG_DISTANCES)

HC=obtainCommTableHC(106)
#HC_f=filteredCommunities(HC,DRUG_DISTANCES)

CL=obtainCommTableCL(0.1925)
#CL_f=filteredCommunities(CL,DRUG_DISTANCES)

heatmapsAlg=function(AP,MCL,PAM,HC,CL,Ref,level,filtered,filename){
  
  if (Ref=="ATCcodes"){
    enrAP=obtainEnrichedSet(AP,"Codes",MergedATCcodes,level)
    enrMCL=obtainEnrichedSet(MCL,"Codes",MergedATCcodes,level)
    enrPAM=obtainEnrichedSet(PAM,"Codes",MergedATCcodes,level)
    enrHC=obtainEnrichedSet(HC,"Codes",MergedATCcodes,level)
    enrCL=obtainEnrichedSet(CL,"Codes",MergedATCcodes,level)
    ATCcodes=Reduce(union, list(enrAP, enrHC, enrPAM,enrMCL,enrCL))
    heatmapcode(ATCcodes,enrAP,enrPAM,enrHC,enrMCL,enrCL,Ref,level,filtered,filename)}
  
  else if (Ref=="MechAction"){
    enrAP_DT=obtainEnrichedSet(AP,"MECHANISM_OF_ACTION",DrugTargets)
    enrMCL_DT=obtainEnrichedSet(MCL,"MECHANISM_OF_ACTION",DrugTargets)
    enrPAM_DT=obtainEnrichedSet(PAM,"MECHANISM_OF_ACTION",DrugTargets)
    enrHC_DT=obtainEnrichedSet(HC,"MECHANISM_OF_ACTION",DrugTargets)
    enrCL_DT=obtainEnrichedSet(CL,"MECHANISM_OF_ACTION",DrugTargets)
    MecAct=Reduce(union, list(enrAP_DT, enrHC_DT, enrPAM_DT,enrMCL_DT,enrCL_DT))
    heatmapcode(MecAct,enrAP_DT,enrPAM_DT,enrHC_DT,enrMCL_DT,enrCL_DT,Ref,level,filtered,filename)
  }
  return(TRUE)
}

heatmapcode=function(enriched,enrAP,enrPAM,enrHC,enrMCL,enrCL,Ref,level,filtered,filename){
  df=data.frame(matrix(0, nrow = 5, ncol = length(enriched)))
  colnames(df)=enriched
  rownames(df)=c("AP","PAM","HC","MCL","CL1")
  alg=c("AP","PAM","HC","MCL","CL1")
  enr=list(AP=enrAP,PAM=enrPAM,HC=enrHC,MCL=enrMCL,CL=enrCL)
  for (i in c(1:5)){
    df[which(rownames(df)==alg[i]),enr[[i]]]=1}
  
  png(paste("../../../DissertationPlots/",filename,".png",sep=""),width = 782, height = 277)
  if (Ref=="ATCcodes"){
    pheatmap(as.matrix(df),color=c("white","grey37"),legend = FALSE,show_colnames = FALSE,fontsize_row=16)
    dev.off()}
  else if (Ref=="MechAction"){
    pheatmap(as.matrix(df),color=c("white","grey37"),legend = FALSE,show_colnames = FALSE,fontsize_row=16)
    dev.off()
  }
}

#Plot heatmaps

heatmapsAlg(AP,MCL,PAM,HC,CL,"ATCcodes",4,TRUE,"heatmapALL106_ATC")

heatmapsAlg(AP,MCL,PAM,HC,CL,"MechAction",4,TRUE,"heatmapALL106_MAct")


#Distribution - Filtered
distribTable=function(comm,alg){
  breaks=c(2,5,10,15,20,25,50,100)
  categories=c("2-4","5-9","10-14","15-19","20-24","25-49","50-99","+100")
  tab=data.frame(algorithm=rep(alg,8),Category=categories,Frequency=numeric(8))
  for (i in c(1:(length(breaks)-1))){
    tab$Frequency[i]=length(which(table(comm$cID)>=breaks[i]& table(comm$cID)<breaks[i+1]))
  }
  tab$Frequency[i+1]=length(which(table(comm$cID)>=breaks[i+1]))
  return(tab)
}

distAP=distribTable(AP,"AP")
distHC=distribTable(HC,"HC")
distMCL=distribTable(MCL,"MCL")
distPAM=distribTable(PAM,"PAM")
distCL=distribTable(CL,"CL1")

dist=rbind(distAP,distHC,distMCL,distPAM,distCL)

categories=c("2-4","5-9","10-14","15-19","20-24","25-49","50-99","+100")

ggplot(dist, aes(x=Category,y=Frequency, fill=algorithm)) +
  geom_bar(stat="identity",position=position_dodge())+colfillScale+theme_Publication()+ scale_x_discrete(limits = categories)+labs(x="Number of Drugs",y="Frequency of clusters",fill="")
ggsave("../../../DissertationPlots/ComparisonDistributionAllFiltered.png")

#Compare 3rd level ATC enrichment
Partit=rbind(HCsevKs_level4,APsevPs_level4,MCLsevIs_level4,PAMsevKs_level4,CLsevDs_level4)
HCsevKs_level4_f$K=HCsevKs_level4$K
MCLsevIs_level4_f$K=MCLsevIs_level4$K
PAMsevKs_level4_f$K=PAMsevKs_level4$K
APsevPs_level4_f$K=APsevPs_level4$K
CLsevDs_level4_f$K=CLsevDs_level4$K
Partit_f=rbind(HCsevKs_level4_f,APsevPs_level4_f,MCLsevIs_level4_f,PAMsevKs_level4_f,CLsevDs_level4_f)

randAP=tidyParamOptimResTable("../DATA/ParameterOptimization/Rand_AP_koptimiz4.csv","AP")
randHC=tidyParamOptimResTable("../DATA/ParameterOptimization/Rand_HC_koptimiz4.csv","HC")
randMCL=tidyParamOptimResTable("../DATA/ParameterOptimization/Rand_MCL_koptimiz4.csv","MCL")
randPAM=tidyParamOptimResTable("../DATA/ParameterOptimization/Rand_PAM_koptimiz4.csv","PAM")
randCL=tidyParamOptimResTable("../DATA/ParameterOptimization/Rand_CL_koptimiz4.csv","CL1")
random=rbind(randAP,randHC,randMCL,randPAM,randCL)
random$random=TRUE

#Remove unnecessary columns from partition tables
Partit=Partit[,c("K","EnrichedClusters","Percentage","EnrichedClustersDT","PercentageDT","Pval_ATC","Pval_DT","algorithm")]
Partit_f=Partit_f[,c("K","EnrichedClusters","Percentage","EnrichedClustersDT","PercentageDT","Pval_ATC","Pval_DT","algorithm")]

Partit$random=FALSE
Partit_f$random=FALSE

#join evertyhing
All3rdLevel=rbind(Partit,random)
All3rdLevel_f=rbind(Partit_f,random)
compareEnrichedClustersRandom=function(AlgComparis,compar,ylab){
  if(compar=="ATC"){
  subset=AlgComparis[,c("K","EnrichedClusters","algorithm","random")]}
  else if (compar == "DT"){
    subset=AlgComparis[,c("K","EnrichedClustersDT","algorithm","random")] 
  }
  subset.m=melt(subset,id=c("algorithm","K","random"))
 p=ggplot(data=subset.m, aes(x=K, y=value,shape=random,colour=algorithm))+
    geom_line(size=1.3,aes(linetype=random, color=algorithm))+labs(x="Number of Clusters",y = ylab)+labs(colour = "",linetype="")+xlim(0,400)+colScale+ geom_vline(xintercept=106,linetype="dashed")  +theme_Publication()+guides(linetype=FALSE)
  return(p)}

highPerc=sapply(c("HC","AP","PAM","CL1","MCL"),function(alg){
  alg_subset=Partit[which(Partit$algorithm==alg & Partit$random==FALSE),]
  return(list(K_ATC=alg_subset$K[which.max(alg_subset$Percentage)],Enr_ATC=alg_subset$EnrichedClusters[which.max(alg_subset$Percentage)],
         K_DT=alg_subset$K[which.max(alg_subset$PercentageDT)],Enr_DT=alg_subset$EnrichedClustersDT[which.max(alg_subset$PercentageDT)]))
})

p1=compareEnrichedClustersRandom(All3rdLevel,"ATC","Enriched Clusters (ATC -3rd level)")+scale_y_continuous(breaks=seq(0,50,10))
p2=compareEnrichedClustersRandom(All3rdLevel,"DT","Enriched Clusters (MoA)")

g=grid.arrange(p1,p2,ncol=2)
ggsave("../../../DissertationPlots/ComparisonALLATC3DTRandomNF.png",g,width=20,height=12,units = "cm")
p3=compareEnrichedClustersRandom(All3rdLevel_f,"EnrichedClusters","Enriched Clusters (ATC -3rd level)")
p4=compareEnrichedClustersRandom(All3rdLevel_f,"EnrichedClustersDT","Enriched Clusters (MoA)")

g=grid.arrange(p3,p4,ncol=2)
ggsave("../../../DissertationPlots/ComparisonALLATC3DTRandomF.png",g,width=20,height=12,units = "cm")

All3rdLevel$Pval_ATC=-log(All3rdLevel$Pval_ATC)
All3rdLevel$Pval_DT=-log(All3rdLevel$Pval_DT)
p3=compareEnrichedClustersRandom(All3rdLevel,"Pval_ATC","-log(pval ATC)")
p4=compareEnrichedClustersRandom(All3rdLevel,"Pval_DT","-log(pval Mech.Action)")

grid.arrange(p3,p4,ncol=2)

dev.off()
