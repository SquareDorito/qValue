library(qvalue)

setwd('/Users/kennoh/Documents/GitHub/qValue/qValue/offline')

timepts=4
is_silac=0
BHfdr=1

p_df = read.table('phospho_df.txt',sep=',',header=TRUE,check.names=FALSE)
link_table = read.table('p_to_np.txt',sep='\t',header=TRUE,check.names=FALSE)
link_df=data.frame(link_table)
rownames(link_df)<-link_df[,2]

stats_df<-data.frame()
fileNames<-character()
pVals<-numeric()
qVals<-numeric()

p_list<-link_df[,2]

for(i in 1:length(p_list)){
  for(j in 1:timepts){
    if(is_silac){
      #silac stuff here
      tempFileName=paste("silac_df/","sample",colnames(p_df)[i],"_timepoint",j,".txt",sep="")
      if(!file.exists(tempFileName)){
        next
      }
    }else{
      #label free stuff
      tempFileName=paste("labelfree_df/","sample",colnames(p_df)[i],"_timepoint",j,".txt",sep="")
      if(!file.exists(tempFileName)){
        next
      }
      dat=read.table(tempFileName,sep='\t',header=TRUE)
      model1 = lm(log(peakArea) ~ channel + np, data = dat, na.action = na.omit)
      model2 = lm(log(peakArea) ~ channel*np, data = dat, na.action = na.omit)
      a = anova(model1, model2)
      
      fileNames<-c(fileNames,tempFileName)
      pVals<-c(pVals,a[['Pr(>F)']][2])
    }
  }
}

if(BHfdr==1){
  print(paste("Controlling FDR with the Benjamin-Hochberg method.\n",sep=""))
  qobj<-qvalue(p=pVals,lambda=0)
}else{
  qobj<-qvalue(p=pVals)
}
qvalues <- qobj$qvalues

if(length(qvalues)==0){
  qobj <- qvalue(pvals, pi0.method="bootstrap")
  qvalues <- qobj$qvalues
}

upper=0.95
while(length(qvalues)==0){
  qobj <- qvalue(pvals, lambda=seq(0, " + upper + ", .01))
  qvalues <- qobj$qvalue
  upper=upper-0.01
  if(upper<0.01){
    break
  }
}

lower=0;
while(length(qvalues)==0){
  qobj <- qvalue(pvals, lambda=seq(" + lower + ", .99, .01))
  qvalues<-qobj$qvalues
  lower=lower+0.01
  if(lower>0.98){
    break
  }
}

if(length(qvalues)==0){
  qobj<-qvalue(pvals,lambda=0)
  qvalues<-qobj$qvalues
}

qVals<-c(qVals,qvalues)

phist<-"phist.png"
png(filename=phist)
hist(pVals,breaks=100)
dev.off()
qplot<-"qplot.png"
png(filename=qplot)
qplot=plot(qobj)
dev.off()
data<-cbind(pVals,qVals)

stats_df <- data.frame(fileNames,pVals,qVals)
write.table(stats_df,file="master_stats_df.txt",sep="\t",row.names=FALSE,quote=FALSE)