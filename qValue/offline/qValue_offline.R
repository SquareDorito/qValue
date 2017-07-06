parseTextFile = function(){
  setwd('C:/Users/knoh1/Documents/qValue_ken/offline')
  link_table = read.table('p_to_np.txt',sep='\t',header=TRUE,check.names=FALSE)
  link_df=data.frame(link_table)
  rownames(link_df)<-link_df[,2]
  write.table(link_df,file="link_df.txt",sep="\t")
  return(link_df)
}

create_lf_df = function(p_df,np_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error){
  link_df<-parseTextFile()
  
  setwd('C:/Users/knoh1/Documents/qValue_ken/offline/labelfree_df')
  
  p_num<-ncol(p_df)
  
  for(i in 2:p_num){
    for (tp in 1:timepts){
      answer <- data.frame()
      id<-character()
      np<-numeric()
      channel<-numeric()
      peakArea<-numeric()
      for (ch in 1:3){
        for (rp in 1:replicates){
          counter=colnames(p_df)[i]
          tempRN = paste("ch", toString(ch), " rep", toString(rp), " timepoint", toString(tp), sep = "")
          if(is.na(p_df[tempRN,counter])){
            next
          }
          tempID = paste("pep",colnames(p_df)[i],"_rep", toString(rp), sep = "")
          channel<-c(channel,ch)
          id<-c(id,tempID)
          np<-c(np,0)
          peakArea<-c(peakArea,p_df[tempRN,counter])
          if(!any(row.names(link_df) == counter)){
            next
          }
          rawNP<-as.character(link_df[counter,3])
          np_list<-as.list(strsplit(rawNP,",")[[1]])
          if(!length(np_list)){
            next
          }
          for(j in 1:length(np_list)){
            if(is.na(np_df[tempRN,np_list[[j]]])){
              next
            }
            channel<-c(channel,ch)
            id<-c(id,tempID)
            np<-c(np,1)
            peakArea<-c(peakArea,np_df[tempRN,np_list[[j]]])
          }
        }
      }
      if(length(peakArea)!=0){
        answer <- data.frame(id,channel,np,peakArea)
        outfile=paste("sample",colnames(p_df)[i],"_timepoint",toString(tp),".txt",sep="")
        write.table(answer,file=outfile,sep="\t",row.names=FALSE,quote=FALSE)
      }
      rm(answer)
      rm(id)
      rm(np)
      rm(channel)
      rm(peakArea)
    }
  }
}

create_silac_df = function(p_df,np_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error){
  link_df<-parseTextFile()
  
  setwd('C:/Users/knoh1/Documents/qValue_ken/offline/silac_df')
  
  p_num<-ncol(p_df)
  
  for(i in 2:p_num){
    for (tp in 1:timepts){
      answer <- data.frame()
      id<-character()
      np<-numeric()
      ratio<-numeric()
      for (rp in 1:replicates){
        counter=colnames(p_df)[i]
        num_rn = paste("ch", toString(1), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        den_rn = paste("ch", toString(2), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        if(is.na(p_df[num_rn,counter]) || is.na(p_df[den_rn,counter])){
          next
        }
        tempID = paste("pep",colnames(p_df)[i],"_rep", toString(rp), sep = "")
        tempRatio=(p_df[num_rn,counter])/(p_df[den_rn,counter])
        id<-c(id,tempID)
        np<-c(np,0)
        ratio<-c(ratio,tempRatio)
        if(!any(row.names(link_df) == counter)){
          next
        }
        rawNP<-as.character(link_df[counter,3])
        np_list<-as.list(strsplit(rawNP,",")[[1]])
        if(!length(np_list)){
          next
        }
        for(j in 1:length(np_list)){
          np_counter<-np_list[[j]]
          cat(paste("i: ",i," np_counter: ",np_counter," peakArea_num: ",np_df[num_rn,np_counter],"\n"))
          if(is.na(np_df[num_rn,np_list[[j]]])){
            next
          }
          if(is.na(np_df[den_rn,np_list[[j]]])){
            next
          }
          id<-c(id,tempID)
          np<-c(np,1)
          np_ratio=(np_df[num_rn,np_list[[j]]])/(np_df[den_rn,np_list[[j]]])
          ratio<-c(ratio,np_ratio)
        }
      }
      if(length(ratio)!=0){
        answer <- data.frame(id,np,ratio)
        outfile=paste("sample",colnames(p_df)[i],"_timepoint",toString(tp),".txt",sep="")
        write.table(answer,file=outfile,sep="\t",row.names=FALSE,quote=FALSE)
      }
      rm(answer)
      rm(id)
      rm(np)
      rm(ratio)
    }
  }
}

library(XML)
## The following are arguments that must be passed
p_filePath = "C:\\Users\\knoh1\\Documents\\qValue_ken\\phospho_data.xml"
np_filePath = "C:\\Users\\knoh1\\Documents\\qValue_ken\\unphospho_data_short.xml"
waitPath = "C:\\Users\\knoh1\\Documents\\qValue_ken\\wait.txt"
replicates = 5
timepts = 12
min_rep = 3
## The following are optional arguments that may be passed
label_free = 2
unpaired = 1
silac_num = 3
silac_den = 2
bh = 0
minPA = 1000
maxPA = Inf
is_min = TRUE
log_LF = FALSE
log_silac = FALSE
is_error = FALSE

## The first two arguments when called must be a path to the XML file and a path to the wait file.
#p_filePath = args[1]
#np_filePath = args[2]
#waitPath = args[2]
## This for loop goes through the rest of the arguments and splits each element into a 1x2 matrix where the elements are
## the sides of the equals sign (ex. arg_mtx = [timepoints, 12]).
## Note that for arguments without an equals sign each element in the matrix will be the same (ex. arg_mtx = [min, min]).


setwd('C:/Users/knoh1/Documents/qValue_ken/offline')
p_df = read.table('phospho_df.txt',sep=',',header=TRUE,check.names=FALSE)
np_df = read.table('unphospho_df.txt',sep=',',header=TRUE,check.names=FALSE)

#if(unpaired==0){
  create_silac_df(p_df,np_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error)
#} else{
  create_lf_df(p_df,np_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error)
#}