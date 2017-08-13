#### Final Plots - Explotatory Data analysis of CMap dataset####

#### Setup ####
library(ConnectivityMap)
library(reshape2)
library(ggplot2)
library("ggsci")
library("ggplot2")
library("gridExtra")
library("grid")
library("ggthemes")
data(instances)
theme_Publication <- function(base_size=14, base_family="helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.key.size= unit(0.6, "cm"),
           legend.margin = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text()
   ))
  
}

#### Retrieving unique instances ####
fact=c("cmap_name","concentration..M.","duration..h.","cell2")
fourFact=instances[,fact]
uniqInstances=unique(fourFact)
nrow(uniqInstances)

#### Plotting Unique instances per drug####
drugs = unique(instances$cmap_name)
retrieveNrInstances=function(drugs,table){
  nrInstancesperDrug=sapply(drugs,function(x){length(which(table$cmap_name==x))})
  nrInstancesperDrugT=table(nrInstancesperDrug)
  morethan7=sum(nrInstancesperDrugT[5:length(nrInstancesperDrugT)])
  nrInstancesperDrugT=nrInstancesperDrugT[1:4]
  nrInstancesperDrugT[5]=morethan7
  names(nrInstancesperDrugT)[5]=">=5"
  return(nrInstancesperDrugT)
}

nrInstancesperDrugT=retrieveNrInstances(drugs,uniqInstances)
df.m=data.frame(Instances=names(nrInstancesperDrugT),value=as.numeric(nrInstancesperDrugT))
df.m$Instances = factor(names(nrInstancesperDrugT),levels=names(nrInstancesperDrugT))
instancesp=ggplot(df.m, aes(Instances, value))+   geom_bar(position = "dodge", stat="identity",color="lightgrey",fill="darkgrey")+xlab("No. Distinct instances")+ylab("No. Drugs")+theme_Publication()+scale_y_continuous(limits=c(0,1500),breaks=seq(0,1500,500))

#### Plotting number of unique instances with distinct cell lines, concentrations and durations per drug#### 
nrperDrug=function(attribute,xlab){
  nrperDrug=sapply(drugs,function(x){length(unique(uniqInstances[which(uniqInstances$cmap_name==x),][,attribute]))})
  nrperDrugT=table(nrperDrug)
  df.m=data.frame(Instances=names(nrperDrugT),value=as.numeric(nrperDrugT))
  instances=ggplot(df.m, aes(Instances, value))+   geom_bar(position = "dodge", stat="identity",color="lightgrey",fill="darkgrey")+xlab(xlab)+theme_Publication()+ylab("")+scale_y_continuous(limits=c(0,1500),breaks=seq(0,1500,500))
  return(instances)
  
}
celllines=nrperDrug("cell2",xlab="No. Cell lines")
duration=nrperDrug("duration..h.",xlab= "No. Treatment durations")
concentration=nrperDrug("concentration..M.", "No. Concentrations")+ylab("No.Drugs")
png("overviewInstancesperDrug.png",width=551,height=355)
grid.arrange(instancesp, celllines, concentration,duration, ncol=2,nrow=2)
dev.off()

#### Plotting number of replicates per instance ####
summaryInst=table(fourFact)
df=data.frame(Replicates=factor(c("0","1","2",">=3"),levels=c("0","1","2",">=3")),Frequency=numeric(4))
for (i in c(1:4)){
  if(i <=3){
  df$Frequency[i]=length(which(summaryInst==i))}
  else{
    df$Frequency[i]=length(which(summaryInst>=i))
  }
}
replicates=ggplot(df, aes(Replicates, Frequency))+   geom_bar(position = "dodge", stat="identity",color="lightgrey",fill="darkgrey")+xlab("Replicates per distinct CMAP instance")+ylab("Frequency")+theme_Publication()

#### Plotting the number of instances per cell line####
cells=table(instances$cell2)
uniqcells=table(uniqInstances$cell2)
df=melt(data.frame("Total"=as.numeric(cells),"Distinct"=as.numeric(uniqcells)))
df$cell=rep(names(cells),2)
#par(mar=c(5.4, 4.5,4.1 ,2.2))
#barplot(cells,horiz=TRUE,col="lightgrey",xlim=c(0,1500),las=2,
#  	main="Cells used for CMAP instances")

instpercell=ggplot(df, aes(cell, value)) +   
geom_bar(aes(fill = variable), position = "dodge", stat="identity")+xlab("Cell lines")+ylab("Number of \nCMap instances")+labs(fill="Instances")+theme_Publication()+theme(legend.position=c(0.8,0.9))
png("Replicates_celllines.png",width=859,height=221)
grid.arrange(replicates,instpercell, ncol=2,nrow=1)
dev.off()
