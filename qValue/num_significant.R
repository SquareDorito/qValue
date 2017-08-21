setwd("C:\\local_qValue")
q_data<-read.table("output.txt",header=FALSE,row.names=1,sep="\t")
counter=0
#cat(paste(nrow(q_data),ncol(q_data)))
for(i in 1:nrow(q_data)){
  for(j in 1:ncol(q_data)){
    #cat(paste(i," ",j,"\n",sep=""))
    if(!is.na(q_data[i,j])){
      if(q_data[i,j]<0.05){
        #cat(q_data[i,j])
        counter=counter+1
      }
    }
  }
}
counter