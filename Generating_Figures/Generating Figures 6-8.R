####R codes for the paper#####
####1. Get the p value of logrank test####
library("survival")

###clinical file downloaded from TCGA and preprocessed for using friendly,our clinical file includes 19 variables.
dat<-read.table("~/your/file/path/clinical1.txt",sep = "\t",header = T,stringsAsFactors = F)

###below is the function used to combine the clinical file and cluster results
fc<-function(clin,cluster){
    colnames(cluster)<-c('Rand','kmean1','kmean2','kmean3','kmean123','spec1','spec2','spec3','spec123','s-spec1','s-spec2','s-spec3','s-spec123','SIM1','SIM2','SIM3','SIM123','Cons1','Cons2','Cons3','Cons123','Cen-M123','Pair-M123','Ker-add123','Our-ed123','Our-edf123','Our-wd123','Our-wde123','Our-wdf123')
    clinccc<-cbind(clin,cluster)
    clinccc$gender<-factor(clinccc$gender)
    clinccc$stat<-as.numeric(clinccc$stat)
    clinccc$time<-as.numeric(clinccc$time)
    clinccc$time<-(clinccc$time)/30
    clinccc$days_to_birth<-as.numeric(clinccc$days_to_birth)
    clinccc$age<-(clinccc$days_to_birth)/365
    for(i in 20:58){
        clinccc[,i]<-factor(clinccc[,i])
    }  ###20:58 indicates the 29 clustering methods
    clinccc
}

###we will get 22 csv files which contained the 29 kinds of clustering results for 22 cancer types
###each file indicates one cancer type

# the order of caner names was corresponding to cancer1.csv to cancer22.csv,which cancer1 means "KIRP"
cancer<-c('BRCA','STAD','LUAD','LUSC','COAD','HNSC','KIRC','BLCA','UVM','PRAD','SARC','KIRP','LIHC','PAAD','ESCA','LGG','MESO','UCEC','THCA','READ','OV','CESC')

for(i in 1:22){
    logrank<-list()
    clin<-dat[which(dat$disease_code==cancer[i]),]
    cluster<-read.table(paste0("~/clustering/cancer",i,".csv"),sep = ",",header = F,stringsAsFactors = F)
    clin$time<-ifelse(clin$time<="0","0.001",clin$time)
    clinccc<-fc(clin,cluster)
    p<-list()
    for(j in 20:58){
        sur<-survdiff(Surv(clinccc$time,clinccc$stat)~clinccc[,j])
        p[[j]]<-1-pchisq(sur$chisq,length(sur$n)-1)
        remove(sur)
        print(p)
    }
    p<-as.data.frame(as.matrix(unlist(p[20:58])))
    colnames(p)<-c("Pvalue")
    p$method<-colnames(clinccc[,c(20:58)])
    p$dis<-rep(cancer[i],29) #9 indicates the 9 clustering methods
    logrank[[i]]<-p ##logrank is the list file which include logrank p values for 22 cancer types
    rm(p)
}

####2. Survival curves and Box plots####
####three examples were provided for ploting the survival curves and box plots

library("survminer")

dat<-read.table("~/your/file/path/clinical2.txt",sep = "\t",header = T,stringsAsFactors = F)

fc2<-function(clin,cluster){
    colnames(cluster)<-c('Rand','kmean1','kmean2','kmean3','kmean123','spec1','spec2','spec3','spec123','s-spec1','s-spec2','s-spec3','s-spec123','SIM1','SIM2','SIM3','SIM123','Cons1','Cons2','Cons3','Cons123','Cen-M123','Pair-M123','Ker-add123','Our-ed123','Our-edf123','MKerWA','Our-wde123','Our-wdf123')
    cluster<-cluster[c("MKerWA")]
    clinccc<-cbind(clin,cluster)
    clinccc$gender<-factor(clinccc$gender)
    clinccc$stat<-as.numeric(clinccc$stat)
    clinccc$time<-as.numeric(clinccc$time)
    clinccc$time<-(clinccc$time)/30
    clinccc$days_to_birth<-as.numeric(clinccc$days_to_birth)
    clinccc$age<-(clinccc$days_to_birth)/365
    clinccc$MKerWA<-factor(clinccc$MKerWA)
    clinccc
}

###plot for Fig6
clin<-dat[which(dat$disease_code=="BRCA"),]
cluster<-read.table("~/your/file/path/cancer1.csv",sep = ",",header = F,stringsAsFactors = F)
clin$time<-ifelse(clin$time<="0","0.001",clin$time)

clinbrca<-fc2(clin,cluster)

fit<-survfit(Surv(time,stat.x)~PAM50,data =clinbrca)
ggsurvplot(fit,data=clinbrca,size=1.5,risk.table=F,font.legend=c(14,"bold","black"),font.y = c(14, "bold", "black"),font.x = c(14, "bold", "black"),ggtheme = theme_minimal(),font.tickslab  = c(14, "bold", "black"),pval = TRUE,conf.int = F,xlab="Time(months)",
legend.title="",legend.labs=c("Basal","Her2","LumA","LumB","Normal"),legend=c(0.89,0.89),palette=c("blue","red","grey","forestgreen","black"),
title=(main="B"),font.title=c(20, "bold", "black"),pval.size=6,pval.coord=c(170,0.75))


fit<-survfit(Surv(time,stat.x)~MKerWA,data =clinbrca)
ggsurvplot(fit,data=clinbrca,size=2,risk.table=F,font.legend=c(12,"bold","black"),font.y = c(12, "bold", "black"),font.x = c(12, "bold", "black"),ggtheme = theme_minimal(),font.tickslab  = c(12, "bold", "black"),pval = TRUE,conf.int = F,xlab="Time(months)",
legend.title="",legend.lab=c("Cluster1","Cluster2","Cluster3","Cluster4"),legend=c(0.89,0.88),palette=c("red","grey","blue","forestgreen"),
title=(main="C"),font.title=c(20, "bold", "black"),pval.size=6,pval.coord=c(170,0.77))


counts_lp<-table(clinbrca$PAM50,clinbrca$MKerWA)
barplot(counts_lp,
xlab="MKerW-A clusters", ylab="Frequency",col=c("darkblue","red","grey","forestgreen","black"),
beside = F,
cex.axis = 1.3,
font.main=2,
cex.main=1.3,
cex.lab=1.3,
cex.names = 1.3,
legend.text=c("Basal","Her2","LumA","LumB","Normal"),
args.legend=list(
x=ncol(counts_lp)+1,
y=max(colSums(counts_lp))+5,
bty = "n", cex=1.3)
)
mtext("A",adj = 0,cex=1.8)



###plot for Fig7

clin<-dat[which(dat$disease_code=="LGG"),]
cluster<-read.table("~/your/file/path/cancer16.csv",sep = ",",header = F,stringsAsFactors = F)
clin$time<-ifelse(clin$time<="0","0.001",clin$time)

clinlgg<-fc2(clin,cluster)


fit<-survfit(Surv(time,stat)~grade,data =clinlgg)
ggsurvplot(fit,data=clinlgg,size=1.5,risk.table=F,font.legend=c(14,"bold","black"),font.y = c(14, "bold", "black"),font.x = c(14, "bold", "black"),ggtheme = theme_minimal(),font.tickslab  = c(14, "bold", "black"),pval = TRUE,conf.int = F,xlab="Time(months)",
legend.title="",legend.labs=c("HIGH","LOW"),legend=c(0.89,0.89),palette=c("forestgreen","black"),
title=(main="B"),font.title=c(20, "bold", "black"),pval.size=6,pval.coord=c(160,0.80))


fit<-survfit(Surv(time,stat)~MKerWA,data =clinlgg)
ggsurvplot(fit,data=clinlgg,size=2,risk.table=F,font.legend=c(12,"bold","black"),font.y = c(12, "bold", "black"),font.x = c(12, "bold", "black"),ggtheme = theme_minimal(),font.tickslab  = c(12, "bold", "black"),pval = TRUE,conf.int = F,xlab="Time(months)",
legend.title="",legend.lab=c("Cluster1","Cluster2","Cluster3","Cluster4"),legend=c(0.89,0.88),palette=c("red","black","blue","forestgreen"),
title=(main="C"),font.title=c(20, "bold", "black"),pval.size=6,pval.coord=c(160,0.75))


counts_lp<-table(clinlgg$grade,clinlgg$MKerWA)
barplot(counts_lp,
xlab="MKerW-A clusters", ylab="Frequency",col=c("forestgreen","black"),
beside = F,
cex.axis = 1.3,
cex.main=1.3,
cex.lab=1.3,
cex.names = 1.3,
legend.text=c("High","Low"),
args.legend=list(
x=ncol(counts_lp)+1,
y=max(colSums(counts_lp)),
bty = "n",cex=1.3))
mtext("A",adj = 0,cex=1.8)


###plot for Fig8
clin<-dat[which(dat$disease_code=="PAAD"),]
cluster<-read.table("~/Documents/MATLAB/cancer15.csv",sep = ",",header = F,stringsAsFactors = F)
clin$time<-ifelse(clin$time<="0","0.001",clin$time)
clinccc<-fc2(clin,cluster)
clinccc<-clinccc[which(clinccc$target_therapy=="YES"|clinccc$target_therapy=="NO"),]

temp<-clinccc[which(clinccc$MKerWA==1),]
fit<-survfit(Surv(time,stat)~target_therapy,data =temp)
ggsurvplot(fit,data=temp,size=1.5,risk.table=F,font.legend=c(16,"bold","black"),font.y = c(16, "bold", "black"),font.x = c(16, "bold", "black"),ggtheme = theme_minimal(),font.tickslab  = c(16, "bold", "black"),conf.int = F,xlab="Time(months)",
legend.title="Cluster1",legend.lab=c("NO","YES"),legend=c(0.85,0.88),
palette=c("forestgreen","black"),pval.size=6,pval=0.68,title=(main="A"),font.title=c(20, "bold", "black"),
pval.coord=c(65,0.8))


temp<-clinccc[which(clinccc$MKerWA==2),]
fit<-survfit(Surv(time,stat)~target_therapy,data =temp)
ggsurvplot(fit,data=temp,size=1.5,risk.table=F,font.legend=c(16,"bold","black"),font.y = c(16, "bold", "black"),font.x = c(16, "bold", "black"),ggtheme = theme_minimal(),font.tickslab  = c(16, "bold", "black"),conf.int = F,xlab="Time(months)",
legend.title="Cluster2",legend.lab=c("NO","YES"),legend=c(0.85,0.88),title=(main="B"),font.title=c(20, "bold", "black"),
palette=c("forestgreen","black"),
pval.size=6,pval=1,pval.coord=c(83,0.8))


temp<-clinccc[which(clinccc$MKerWA==3),]
fit<-survfit(Surv(time,stat)~target_therapy,data =temp)
ggsurvplot(fit,data=temp,size=1.5,risk.table=F,font.legend=c(16,"bold","black"),font.y = c(16, "bold", "black"),font.x = c(16, "bold", "black"),ggtheme = theme_minimal(),font.tickslab  = c(16, "bold", "black"),conf.int = F,xlab="Time(months)",
legend.title="Cluster3",legend.lab=c("NO","YES"),legend=c(0.85,0.88),title=(main="C"),font.title=c(20, "bold", "black"),
palette=c("forestgreen","black"),
pval.size=6,pval="p < 0.001",pval.coord=c(65,0.8))


temp<-clinccc[which(clinccc$MKerWA==4),]
fit<-survfit(Surv(time,stat)~target_therapy,data =temp)
ggsurvplot(fit,data=temp,size=1.5,risk.table=F,font.legend=c(16,"bold","black"),font.y = c(16, "bold", "black"),font.x = c(16, "bold", "black"),ggtheme = theme_minimal(),font.tickslab  = c(16, "bold", "black"),conf.int = F,xlab="Time(months)",
legend.title="Cluster4",legend.lab=c("NO","YES"),legend=c(0.85,0.88),title=(main="D"),font.title=c(20, "bold", "black"),
palette=c("forestgreen","black"),
pval.size=6,pval="p < 0.001",pval.coord=c(49,0.8))
