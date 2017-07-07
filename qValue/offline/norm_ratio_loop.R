library(qvalue)

setwd('C:/Users/knoh1/Documents/qValue_ken/offline')

timepts=4
is_silac=0

p_df = read.table('phospho_df.txt',sep=',',header=TRUE,check.names=FALSE)

stats_df<-data.frame()
fileNames<-character()
pVals<-numeric()
qVals<-numeric()

for(i in 2:ncol(p_df)){
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

qobj<-qvalue(p=pVals)
qvalues <- qobj$qvalues
qVals<-c(qVals,qvalues)
stats_df <- data.frame(fileNames,pVals,qVals)
write.table(stats_df,file="master_stats_df.txt",sep="\t",row.names=FALSE,quote=FALSE)