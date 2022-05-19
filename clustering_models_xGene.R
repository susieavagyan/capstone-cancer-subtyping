###xGene
library(dplyr)
library(readr)
library(GO.db)
library(org.Hs.eg.db)
library("GOSemSim")
library("jqr")
library(dplyr)
library(tidyr)


######### Calculate distance matrix between tumor samples

dataName="./data/mut_cnv_onehot.csv"
dat=read.csv(dataName, header=T)


FSM=read.table("results/xgene/FSM.txt", header=F, sep="\t", skip=1)
FSM[which(FSM==NA)]=0

WM=read.table("results/xgene/WM.txt", header=F, sep="\t", skip=1)

dat1 <- dat %>% select_if(colnames(dat) %in% WM$V1)
dat2 <-  dat %>% dplyr::select(patient_id)

dat <- cbind(dat2, dat1)

WM = WM[WM$V1 %in% colnames(dat),]

FSM = FSM[FSM$V1 %in% colnames(dat),]

p=nrow(dat)
e=ncol(dat)
w=nrow(WM)


m=match(gsub("-",".",as.character(WM[,1])),names(dat), nomatch=0)
SMIM=matrix(as.integer(dat[,m]!="0"),p,w)


M2=WM[,-1]*FSM[,-1]
M2=as.matrix(M2)
M2=as.vector(M2)
ix=which(is.na(M2)==T)
M2[ix]=0
M2=matrix(M2,w,w)

ix=which(colSums(SMIM)!=0)
SMIM=SMIM[,ix]
M2=M2[ix,ix]

D=matrix(0,p,p)
ng=2                           #the number of driver genes to be considered

for(i in 1:(p-1)){ #print (i)
  for(j in (i+1):p){
    
    r=which(SMIM[i,]==1)
    c=which(SMIM[j,]==1)
    lc=length(c)
    lr=length(r)
    
    s=rep(0,ng)
    
    if(min(length(c),length(r))>0){
      
      M3=matrix(0,lr+5,lc+5)
      M3[1:lr,1:lc]=M2[r,c]
      
      for(l in 1:ng){
        s[l]=max(M3,na.rm=T);ix=which(M3==s[l],arr.ind=T)
        M3=M3[-ix[1],-ix[2]]
      } 
      
    }  
    
    D[i,j]=sum(s)/ng
  }
}

D=1-(D+t(D))
diag(D)=0
#plot(hclust(as.dist(D),method="ward.D"))

row.names(D)=as.character(dat$bcr_patient_barcode)

outputName="results/xgene/distance.txt"
write.table(D, file=outputName, col.names=F, row.names=T, sep="\t", quote=F)


##clustering and evaluation


library(survival)
library(logistf)
clin_patients = read.csv('./data/clinical_patients.csv') ##clinical data
clin_samples = read.csv('./data/clinical_primary_raw.csv')
clin = merge(clin_patients, clin_samples, by = "PATIENT_ID")

dat1 = merge(clin, dat, by.x = "SAMPLE_ID", by.y = "patient_id")

cancerName="Colorectal"
save_dir <- file.path("./results/xgene")
if (!dir.exists(save_dir)){
  dir.create(save_dir, recursive = TRUE)
} 

distFileName="./results/xgene/distance.txt"

figureName="./results/xgene/cluster_survival_association.pdf"
clusterLabelName="./results/xgene/cluster_labels.csv"       

overview=as.data.frame(matrix(NA,1,3))
overview[1,1]=cancerName  

nc=5  #the number of clusters                                                                
pdf("./results/xgene/clustering_and_evaluation.pdf",width=9.2, height=6.5)

#race significance
# dat=read.csv(clnFileName, header=T)
dm=read.csv(distFileName, header=F, sep="\t")
dm=dm[,-1]   

race=as.character(dat1$RACE)

my.clust=hclust(as.dist(dm),method="ward.D")
cluster=cutree(my.clust,nc)

write.csv(cbind(dat[1],cluster),file=clusterLabelName,row.names=F,quote=F)

r=which(race=="Asian-far east/indian subcont"|race=="Black or african american"|race=="White")
tb=table(race[r],cluster[r])
row.names(tb)=paste(c("AN","BL","WH")," (",table(race[r]),")",sep="")
print(tb)
p.value=fisher.test(tb[,],simulate.p.value=T)$p.value
print(fisher.test(tb[,],simulate.p.value=T))


## survival significance 

p.contr=1
contr="--"

for(c1 in 1:nc){
  for(c2 in c1:nc){
    my.surv.contr=coxph(Surv(OS_MONTHS,OS_STATUS=="1:DECEASED")~AGE_AT_SEQUENCING + as.factor(as.integer(x %in% c(c1,c2))), data=dat1)
    if(summary(my.surv.contr)$coefficient[2,5]<p.contr){
      p.contr=summary(my.surv.contr)$coefficient[2,5]
      contr=paste("(C",c1," vs ",sep="")
      if(c1!=c2)contr=paste("(C",c1,",",c2," vs ",sep="")
      contr=paste(contr,"others)",sep="")       
    } 
  }
}

x=cluster
my.surv.anova=coxph(Surv(OS_MONTHS,OS_STATUS=="1:DECEASED")~AGE_AT_SEQUENCING + as.factor(as.integer(x)), data=dat1)
p.anova=anova(my.surv.anova)[3,4]

print(summary(my.surv.anova))

xlab="Overall survival months"
ylab="Survival probability"
my.fit=survfit(Surv(OS_MONTHS,OS_STATUS=="1:DECEASED")~as.factor(x),data=dat1)
plot(my.fit, col=colorCodes,lty=5,lwd=1.2, xlab="",ylab=ylab,main="Cluster-survival association",cex=.8, cex.axis=.8)    #main=mf[k,1]
mtext(side=1,xlab, cex=.8, line=2) 

mtext(side=3,paste("Cox-PH p-value:",format(round(min(p.anova),4),scientific=F),contr),cex=.8, at=65,line=.3)
mtext(side=3,paste(cancerName,"- results of xGeneModel", sep=" "),line=-1.5, outer=T, cex=1.2)

list(p.value=format(round(min(p.contr),4),scientific=F))

close.screen(all.screens = T)
split.screen(c(2,1))
split.screen(c(1,3), screen = 1)
split.screen(c(1,nc), screen = 2)
library(dendextend)
colorCodes <- c("red","green","blue","purple","black","orange4")
dend=as.dendrogram(my.clust, hang=4)
labels_colors(dend) <- colorCodes[cluster][order.dendrogram(dend)]
labels(dend)=rep("E",nrow(dat))

height=sort(my.clust$height,decreasing=T)
h=(height[nc]+height[nc-1])/2

#### 
lab1=labels_colors(dend)
#lab2=lab1[1:length(lab1)]
lab2=lab1[2:length(lab1)]
lab2=c(lab1[2:length(lab1)],"A")
mark=lab1==lab2
mark=which(mark=="FALSE")
pos=(mark+c(0,mark[1:(nc-1)]))/2; pos[2]=pos[2]-10
size=mark-c(0,mark[1:(nc-1)])
####

screen(3)
par(mar=c(2.3, 4, 4, 2) + 0.1)
plot(dend, main="Mutation-based clustering",cex.axis=.8)
abline(h=h,lty=2)
mtext(side=1,at=pos,size,line=0.5, cex=0.7)

screen(4)
par(mar=c(4, 4, 4, 2) + 0.1)   #par(mar=c(4, 4, 3, 2) + 0.1)
x=cluster
overview[1,2]=format(round(min(p.contr),4),scientific=F)

screen(5)
par(mar=c(4.5, 4, 4, 2) + 0.1)
barplot(t(tb/rowSums(tb)*100),col=colorCodes,ylab="Percentage partition", main="Cluster-race association",cex.names=.7,cex=.8)
mtext(side=1, paste("Fisher's test p-value:", format(round(p.value,4),scientific=F)),line=2.5,cex=.8)
overview[1,3]=format(round(p.value,4),scientific=F)

######
p=nrow(dat1)
e=ncol(dat1)
b=e-ng+1

for(r in 1:nc){
  set=dat1[which(cluster==r),52:e]!="0"
  ls=sort(colSums(set), decreasing=T)/length(which(cluster==r))
  
  ix=which(ls==0)
  names(ls)[ix]="----"   
  
  screen(r+5)
  par(mar=c(4.5, 4, 3.5, 2) + 0.1)
  barplot(ls[1:10],horiz=T, xlab="Frequency",cex.lab=.9, cex.axis=.8, las=1, 
          cex.names=.6, main=list(paste("cluster",r),col=colorCodes[r]),xlim=c(0,1.0),cex=.8)
  
  
}

dev.off()   



names(overview)=c("cancer","survival_significance", "race_significance")
write.table(overview, file="./results/xgene/result_overview.txt", col.names=T, row.names=F, sep="\t", quote=F)





##calculating silhouette score

library(cluster)
cluster_assignments <- cbind(dat[1],cluster)
assignment_vector <- cluster_assignments$cluster
names(assignment_vector) <- cluster_assignments$patient_id
si <- silhouette(assignment_vector, D )
summary(si)
